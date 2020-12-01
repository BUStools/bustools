#include <iostream>
#include <fstream>
#include <cstring>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_predict.h"
#include <iostream>
#include <fstream>
#include <cmath>

#include <Eigen/Core>
#include <unsupported/Eigen/SpecialFunctions>
#include <LBFGSB.h>
#include <cppoptlib/meta.h>
#include <cppoptlib/boundedproblem.h>
#include <cppoptlib/solver/lbfgsbsolver.h>
#include <thread>
#include <memory>
#include <mutex>

//based on the boost implementation
//r is the size parameter
inline double DensityNegBin(double k, double r, double mu) {
	//so, the standard parameterization is r and p, but we use mean (mu) and size
	//mean = mu = r(1-p)/p <=> p = r/(mu+r), size (i.e. r) in this parameterization and r in the standard one is the same.
	double p = r / (mu + r);
	return exp(lgamma(r + k) - lgamma(r) - lgamma(k + 1)) * pow(p, r) * pow((1 - p), k);
}


//So, this one is different from the negative binomial below in that the zero is not included.
//We here want to look at the log likelihood for observed molecules (i.e. not including the zeros) and how they would fit the negative
//binomial, which means that the density function needs to be modified. The probability for zero needs to be set to zero, and the other
//values scaled up to get to a sum of 1.
inline double ZTNBLogLikelihood(const double* hist, size_t histLen, double size, double mu) {
	//first, get the probability of zero, probZero - all other values should be divided by 1-probZero to normalize the sum of non-zeros to 1
	double probZero = DensityNegBin(0, size, mu);
	double scale = 1 - probZero;

	//now calculate the log likelihood
	double ll = 0;
	for (size_t i = 0; i < histLen; ++i) {
		ll += log(DensityNegBin(double(i + 1), size, mu) / scale) * (*(hist + i));
	}

	return ll;
}

inline double NegBinLogLikelihood(const double* hist, size_t histLen, double size, double mu, double zeroCount) {
	//loglikelihood for nonzero counts
	double ll = 0;
	for (size_t i = 0; i < histLen; ++i) {
		ll += log(DensityNegBin(double(i + 1), size, mu)) * (*(hist + i));
	}
	//...and for the zero counts
	ll += log(DensityNegBin(0, size, mu)) * zeroCount;
	
	return ll;
}

//using LBFGSpp
class Optimizer
{
public:
	//make the variables public to be able to update them without a lot of extra code
	double mean = 0;
	double zeroCounts = 0;
	double totMol = 0;
	const double* hist = nullptr;
	size_t histLen = 0;
	Optimizer() {}

	//returns objective function value. The grad vector should be modified (i.e. set gradients)
	//the vector x is the in parameters, in this case the size (i.e. r) parameter in the distribution
	//The math here is the same as in preseqR
	double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
		double firstTerm = Eigen::numext::digamma(x[0]) * zeroCounts;
		for (size_t i = 0; i < histLen; ++i) {
			firstTerm += *(hist+i) * Eigen::numext::digamma(double(i + 1 + x[0]));
		}
		firstTerm /= totMol;
		
		double secondTerm = Eigen::numext::digamma(x[0]);

		double thirdTerm = log(x[0]) - log(x[0] + mean);

		grad[0] = -firstTerm + secondTerm - thirdTerm;

		return -NegBinLogLikelihood(hist, histLen, x[0], mean, zeroCounts) / totMol;
	}
};

//using CppOptimizationLibrary
class Optimizer2 : public cppoptlib::BoundedProblem<double, 1> {
public:
	//make the variables public to be able to update them without a lot of extra code
	double mean = 0;
	double zeroCounts = 0;
	double totMol = 0;
	const double* hist = nullptr;
	size_t histLen = 0;

	//The math here is the same as in preseqR

	double value(const TVector& x) {
		return -NegBinLogLikelihood(hist, histLen, x[0], mean, zeroCounts) / totMol;
	}

	void gradient(const TVector& x, TVector& grad)
	{
		double firstTerm = Eigen::numext::digamma(x[0]) * zeroCounts;
		for (size_t i = 0; i < histLen; ++i) {
			firstTerm += *(hist + i) * Eigen::numext::digamma(double(i + 1 + x[0]));
		}
		firstTerm /= totMol;

		double secondTerm = Eigen::numext::digamma(x[0]);

		double thirdTerm = log(x[0]) - log(x[0] + mean);

		grad[0] = -firstTerm + secondTerm - thirdTerm;
	}
};

//size and mu are both in and out parameters - the in parameters are the inital values in the expectation maximization algorithm
//the hist is just a pointer into the large hist vector to avoid unnecessary copying
//returns the log likelihood
//This function uses LBFGSpp
double PredictZTNBEmAlg1(const double* hist, size_t histLen, double& size, double& mu) {

	double zeroProbability = DensityNegBin(0, size, mu);
	double histSum = 0; //S in the R code
	for (size_t i = 0; i < histLen; ++i) {
		histSum += *(hist + i);
	}
	//estimate the total number of molecules
	double totMol = histSum / (1 - zeroProbability); //L in the R code

	//the estimated number of zero counts
	double zeroCounts = totMol * zeroProbability;

	//estimate the mean and variance
	double mean = 0, variance = 0;
	double precalcMeanSum = 0; //reuse this sum in the EM loop below since it is constant throughout the algorithm
	for (size_t i = 0; i < histLen; ++i) {
		precalcMeanSum += *(hist + i) * double(i + 1);//remember that the hist index is off by one, i.e. that the hist[0] corresponds to molecules with 1 copy and so forth
	}
	mean = precalcMeanSum / totMol;
	//..and the variance
	for (size_t i = 0; i < histLen; ++i) {
		variance += *(hist + i) * pow(double(i + 1) - mean, 2);
	}
	variance = (variance + mean * mean * zeroCounts) / (totMol - 1);

	Optimizer op;

	// set all params
	op.mean = mean;
	op.zeroCounts = zeroCounts;
	op.totMol = totMol;
	op.hist = hist;
	op.histLen = histLen;


	//using LBFGSpp
	LBFGSpp::LBFGSBParam<double> param;
	//So, the epsilon and max iterations here are for the LBFGS algorithm, not the EM, that comes later!
	param.epsilon = 1e-8; //It seems that the optim function in R uses about 1e-8
	param.max_iterations = 100; //defaults to 100 in the optim function in R
	param.max_linesearch = 500; //don't want it to fail...

	// Create solver and function object
	LBFGSpp::LBFGSBSolver<double> solver(param);

	// Initial guess
	Eigen::VectorXd x = Eigen::VectorXd::Zero(1);

	if (variance > mean) {
		x[0] = mean * mean / (variance - mean);
	}
	else {
		x[0] = size;
	}

	// Bounds
	Eigen::VectorXd lb = Eigen::VectorXd::Constant(1, 0.0001);
	Eigen::VectorXd ub = Eigen::VectorXd::Constant(1, 10000);

	// x will be overwritten to be the best point found
	double fx = 0;
	int numIter = solver.minimize(op, x, fx, lb, ub);


	size_t iter = 0;
	double lastNegLL = 10000000000000;
	double currNegLL = -ZTNBLogLikelihood(hist, histLen, x[0], mean);


	//termination criteria for the EM algorithm
	const double MAX_ERROR = 1e-8;
	const double MAX_ERROR_FAST = 1e-5;
	const size_t MAX_ITER = 100000;
	const size_t ITER_FAST_LIMIT = 200; //if the number of iterations goes over this number, start using the lower error threshold, MAX_ERROR_FAST, to quicken things up


	//The EM algorithm starts here
	while (fabs(lastNegLL - currNegLL) / histSum > MAX_ERROR&&
		iter < MAX_ITER &&
		!(iter >= ITER_FAST_LIMIT && fabs(lastNegLL - currNegLL) / histSum <= MAX_ERROR_FAST))
	{
		lastNegLL = currNegLL;

		//update distribution params
		size = x[0];
		mu = mean;

		// E step: estimate the number of unobserved species

		zeroProbability = DensityNegBin(0, size, mu);
		totMol = histSum / (1 - zeroProbability);
		zeroCounts = totMol * zeroProbability;

		//mean and variance
		mean = precalcMeanSum / totMol;
		//...and variance
		for (size_t i = 0; i < histLen; ++i) {
			variance += *(hist + i) * pow(double(i + 1) - mean, 2);
		}
		variance = (variance + mean * mean * zeroCounts) / (totMol - 1);



		// M step: estimate the parameters size and mu
		//rerun the LBFGSB

		// set params
		op.mean = mean;
		op.zeroCounts = zeroCounts;
		op.totMol = totMol;

		if (variance > mean) {
			x[0] = mean * mean / (variance - mean);
		}
		else {
			x[0] = size;
		}

		//avoid printing of warning text in console
		if (x[0] > 10000.0) {
			x[0] = 10000.0;
		}
		if (x[0] < 0.0001) {
			x[0] = 0.0001;
		}


		int numIter = solver.minimize(op, x, fx, lb, ub);

		currNegLL = -ZTNBLogLikelihood(hist, histLen, x[0], mean);

		//std::cout << "error: " << (lastNegLL - currNegLL) / histSum << "\n";

		++iter;

		//if (iter > 300) break;//tmp
	}

	//std::cout << "iterations: " << iter << " error: " << (lastNegLL - currNegLL) / histSum << "\n";

	return -currNegLL;
}

//Similar to above, but uses CppOptimizationLibrary
double PredictZTNBEmAlg2(const double* hist, size_t histLen, double& size, double& mu) {

	double zeroProbability = DensityNegBin(0, size, mu);
	double histSum = 0; //S in the R code
	for (size_t i = 0; i < histLen; ++i) {
		histSum += *(hist + i);
	}
	//estimate the total number of molecules
	double totMol = histSum / (1 - zeroProbability); //L in the R code

	//the estimated number of zero counts
	double zeroCounts = totMol * zeroProbability;

	//estimate the mean and variance
	double mean = 0, variance = 0;
	double precalcMeanSum = 0; //reuse this sum in the EM loop below since it is constant throughout the algorithm
	for (size_t i = 0; i < histLen; ++i) {
		precalcMeanSum += *(hist + i) * double(i + 1);//remember that the hist index is off by one, i.e. that the hist[0] corresponds to molecules with 1 copy and so forth
	}
	mean = precalcMeanSum / totMol;
	//..and the variance
	for (size_t i = 0; i < histLen; ++i) {
		variance += *(hist + i) * pow(double(i + 1) - mean, 2);
	}
	variance = (variance + mean * mean * zeroCounts) / (totMol - 1);

	Optimizer2 op;
	// set all params
	op.mean = mean;
	op.zeroCounts = zeroCounts;
	op.totMol = totMol;
	op.hist = hist;
	op.histLen = histLen;


	//using CppOptimizationLibrary
	op.setLowerBound(Optimizer2::TVector::Ones(1) * 0.0001);
	op.setUpperBound(Optimizer2::TVector::Ones(1) * 10000);
	cppoptlib::LbfgsbSolver<Optimizer2> solver;
	Optimizer2::TVector x = Optimizer2::TVector::Zero(1);

	if (variance > mean) {
		x[0] = mean * mean / (variance - mean);
	}
	else {
		x[0] = size;
	}

//	solver.setDebug(cppoptlib::DebugLevel::Low);
	// x will be overwritten to be the best point found
	solver.minimize(op, x);


	// Bounds
	Eigen::VectorXd lb = Eigen::VectorXd::Constant(1, 0.0001);
	Eigen::VectorXd ub = Eigen::VectorXd::Constant(1, 10000);



	size_t iter = 0;
	double lastNegLL = 10000000000000;
	double currNegLL = -ZTNBLogLikelihood(hist, histLen, x[0], mean);


	//termination criteria for the EM algorithm
	const double MAX_ERROR = 1e-8;
	const double MAX_ERROR_FAST = 1e-5;
	const size_t MAX_ITER = 100000;
	const size_t ITER_FAST_LIMIT = 200; //if the number of iterations goes over this number, start using the lower error threshold, MAX_ERROR_FAST, to quicken things up


	//The EM algorithm starts here
	while (fabs(lastNegLL - currNegLL) / histSum > MAX_ERROR&&
		iter < MAX_ITER &&
		!(iter >= ITER_FAST_LIMIT && fabs(lastNegLL - currNegLL) / histSum <= MAX_ERROR_FAST))
	{
		lastNegLL = currNegLL;

		//update distribution params
		size = x[0];
		mu = mean;

		// E step: estimate the number of unobserved species

		zeroProbability = DensityNegBin(0, size, mu);
		totMol = histSum / (1 - zeroProbability);
		zeroCounts = totMol * zeroProbability;

		//mean and variance
		mean = precalcMeanSum / totMol;
		//...and variance
		for (size_t i = 0; i < histLen; ++i) {
			variance += *(hist + i) * pow(double(i + 1) - mean, 2);
		}
		variance = (variance + mean * mean * zeroCounts) / (totMol - 1);



		// M step: estimate the parameters size and mu
		//rerun the LBFGSB

		// set params
		op.mean = mean;
		op.zeroCounts = zeroCounts;
		op.totMol = totMol;

		if (variance > mean) {
			x[0] = mean * mean / (variance - mean);
		}
		else {
			x[0] = size;
		}

		//avoid printing of warning text in console
		if (x[0] > 10000.0) {
			x[0] = 10000.0;
		}
		if (x[0] < 0.0001) {
			x[0] = 0.0001;
		}

		solver.minimize(op, x);
		//std::cout << "val: " << x[0] << "\n";

		currNegLL = -ZTNBLogLikelihood(hist, histLen, x[0], mean);

		//std::cout << "error: " << (lastNegLL - currNegLL) / histSum << "\n";

		++iter;

		//if (iter > 300) break;//tmp
	}

	//std::cout << "iterations: " << iter << " error: " << (lastNegLL - currNegLL) / histSum << "\n";

	return -currNegLL;
}

//size and mu are out parameters describing the negative binomial
double PredictZTNBForGene(const double* hist, size_t histLen, double t, double& size, double& mu) {
	double histSum = 0; //S in the R code
	for (size_t i = 0; i < histLen; ++i) {
		histSum += *(hist + i);
	}
	if (histSum == 0.0) {
		return 0.0;//nothing to do...
	}
	//initial values of negative binomial, taken from R code
	size = 1.0;
	mu = 0.5;
	//fit the ZTNB (will update size and mu)
	//So, the trick here is to first use Alg1 - it is faster, but fails sometimes. If it fails,
	//use Alg2
	try {
		PredictZTNBEmAlg1(hist, histLen, size, mu);
	}
	catch (std::exception&)
	{
		PredictZTNBEmAlg2(hist, histLen, size, mu);
	}
	
	//std::cout << "Mu: " << mu << " Size: " << size << "\n";

	//estimate the total number of molecules
	double zeroProbability = DensityNegBin(0, size, mu);
	double totMol = histSum / (1 - zeroProbability); //L in the R code

	//The prediction is based on the assumption that the size parameter remains the same, while
	//the mean is scaled up with t. We then look at how many non-zero molecules we would get
	double zeroProbScaled = DensityNegBin(0, size, mu*t);
	double result = totMol * (1 - zeroProbScaled);
	if (result < histSum) { //safety check, we should never get fewer molecules than what we started with!
		result = histSum;
	}

	return result;
}

//This is a thread scheduler
class PredictionExecuter
{
public:
	PredictionExecuter(const std::vector<double>& hists, const std::vector<size_t>& histLengths, double t, std::vector<double>& predVals, std::vector<double>& sizeVals, std::vector<double>& muVals, const uint32_t histmax) 
		: m_hists(hists)
		, m_histLengths(histLengths)
		, m_t(t)
		, m_predVals(predVals)
		, m_sizeVals(sizeVals)
		, m_muVals(muVals)
		, m_histmax(histmax)
	{}
	void Execute(int numThreads) {
		//create the threads
		for (int i = 0; i < numThreads; ++i) {
			m_threads.push_back(std::shared_ptr<std::thread> (new std::thread(&PredictionExecuter::ThreadFunc, this)));
		}
		//wait for them to finish
		for (int i = 0; i < numThreads; ++i) {
			m_threads[i]->join();
		}
		std::cout << "\n";
	}
private:
	void ThreadFunc() {
		size_t i = 0;
		while (GetNextIndex(i)) {
			double size = 0;
			double mu = 0;
			m_predVals[i] = PredictZTNBForGene(&m_hists[i * m_histmax], m_histLengths[i], m_t, size, mu);
			m_sizeVals[i] = size;
			m_muVals[i] = mu;
		}
	}
	bool GetNextIndex(size_t& index) {
		std::lock_guard<std::mutex> lg(m_mutex);
		if (m_nextIndex == m_histLengths.size()) {
			return false;
		} else {
			std::cout << "\rProcessing gene: " << m_nextIndex + 1 << " of " << m_histLengths.size();//so, we use indexes starting at 1 in the printout
			index = m_nextIndex++;
			return true;
		}
	}
	std::mutex m_mutex;
	size_t m_nextIndex = 0;
	const std::vector<double>& m_hists;
	const std::vector<size_t>& m_histLengths;
	double m_t;
	std::vector<double>& m_predVals;
	std::vector<double>& m_sizeVals;
	std::vector<double>& m_muVals;
	const uint32_t m_histmax;
	std::vector<std::shared_ptr<std::thread>> m_threads;
};

void bustools_predict(Bustools_opt &opt) {

	//Prepare and load histograms
	//////////////////

	std::cout << "Test 1\n";

	//read the hist file
	std::string hist_ifn = opt.predict_input + ".hist.txt";
	std::string counts_ifn = opt.predict_input + ".mtx";
	std::string gene_ifn = opt.predict_input + ".genes.txt";
	std::string barcode_ifn = opt.predict_input + ".barcodes.txt";

	std::string corr_counts_ofn = opt.output + ".mtx";
	std::string nb_params_ofn = opt.output + "nb_params.txt";
	std::string corr_gene_ofn = opt.output + ".genes.txt";
	std::string corr_barcode_ofn = opt.output + ".barcodes.txt";

	//just copy the genes and barcodes files, they don't change
	copy_file(gene_ifn, corr_gene_ofn);
	copy_file(barcode_ifn, corr_barcode_ofn);

	//Get gene list associated with count matrix and histograms
	std::vector<std::string> genes;
	parseGenesList(gene_ifn, genes);

	//Allocate histograms
	//Indexed as gene*histmax + histIndex
	size_t n_genes = genes.size();
	const uint32_t histmax = 100;//set molecules with more than histmax copies to histmax 
	std::vector<double> histograms = std::vector<double>(n_genes * histmax, 0);
	std::vector<size_t> histogramLengths = std::vector<size_t>(n_genes, 0);
	std::string line;

	std::cout << "Test 2\n";

	//load histograms file
	{
		std::ifstream inf(hist_ifn);
		size_t geneIndex = 0;
		while (std::getline(inf, line)) {
			std::stringstream ss(line);
			double num = 0;
			size_t histIndex = 0;
			while (ss >> num) {
				histograms[geneIndex * histmax + histIndex++] = num;
			}
			histogramLengths[geneIndex] = histIndex;
			++geneIndex;
		}
	}

	std::cout << "Test 3\n";

	//fix the histograms if they have only a single non-zero value (i.e. for example looks like this: 1 0 0 0)
	//however, do not fix completely empty histograms
	//the strategy is to always add a one after the last item in the histogram 
	//(this will be next to the only non-zero number, which will be last)
	//the strategy may lead to that a gene gets more counts. This, however, doesn't matter since we only calculate 
	//a scaling factor for the gene, so that will be normalized out.
	for (size_t i = 0; i < histogramLengths.size(); ++i) {
		size_t nonZeros = 0;
		for (size_t j = 0; j < histogramLengths[i]; ++j) {
			if (histograms[i * histmax + j] != 0) {
				++nonZeros;
				if (nonZeros > 1) {
					break;
				}
			}
		}
		if (nonZeros == 1) {
			if (histogramLengths[i] == histmax) { //rare case, add a count before instead of after
				histograms[i * histmax + histmax - 2] = 1;
			} else {
				histograms[i * histmax + histogramLengths[i]] = 1;
				histogramLengths[i]++;//also extend the histogram one step;
			}
			//so, if we can, remove one count from the current value
			if (histograms[i * histmax + histogramLengths[i] - 1] > 1) {
				histograms[i * histmax + histogramLengths[i] - 1] -= 1.0;
			}
		}
	}

	std::cout << "Test 4\n";

	//Predict
	//////////////////
	
	std::vector<double> predVals(n_genes, 0);
	std::vector<double> sizeVals(n_genes, 0);
	std::vector<double> muVals(n_genes, 0);
	int numThreads = std::thread::hardware_concurrency();
	std::cout << "Using " << numThreads << " threads\n";
	PredictionExecuter pe(histograms, histogramLengths, opt.predict_t, predVals, sizeVals, muVals, histmax);
	pe.Execute(numThreads);
	
	//calculate the sum of all histograms and all predvals:
	double histSum = 0;
	std::vector<double> umisPerGene(n_genes, 0);
	std::vector<double> countsPerGene(n_genes, 0);
	double predSum = 0;

	for (size_t i = 0; i < predVals.size(); ++i) {
		for (size_t j = 0; j < histogramLengths[i]; ++j) {
			umisPerGene[i] += histograms[i * histmax + j];
			countsPerGene[i] += histograms[i * histmax + j]*(j+1);
		}
		histSum += umisPerGene[i];
		predSum += predVals[i];
	}

	//so the genes should be scaled according to the following:
	//geneScaling = histSum/predSum * predVal/umisPerGene;
	std::vector<double> geneScaling(n_genes, 1.0);//1.0 for all empty genes
	double globScale = histSum / predSum;
	for (size_t i = 0; i < predVals.size(); ++i) {
		if (umisPerGene[i] > 0) {
			geneScaling[i] = globScale * predVals[i] / umisPerGene[i];
		}
	}
	
	//Modify counts
	//////////////////


	//read and write the counts matrix
	{
		std::ifstream inf(counts_ifn);
		std::string test;
		std::ofstream of(corr_counts_ofn);
		//first number is cells, second genes
		//first handle header, we just leave that untouched
		while (std::getline(inf, line)) {
			of << line << '\n';
			if ((!line.empty()) && line[0] != '%') {
				break;
			}
		}
		//now all data with scaling
		size_t cell = 0, gene = 0;
		double count = 0;
		while (std::getline(inf, line)) {
			if ((!line.empty()) && line[0] == '%') {
				of << line << '\n';
				continue;
			}
			std::stringstream ss(line);
			if (ss >> cell >> gene >> count) {
				of << cell << " " << gene << " " << count * geneScaling[gene-1] << '\n';
			}
		}
	}
	
	//write negative binomial params (file with header)
	{
		std::ofstream of(nb_params_ofn);
		
		//header
		of << "gene\tmu\tsize\tUMIs\tcounts\n";

		for (size_t i = 0; i < genes.size(); ++i) {
			of << genes[i] << '\t' << muVals[i] << '\t' << sizeVals[i] << '\t' << umisPerGene[i] << '\t' << countsPerGene[i] << '\n';
		}
	}
	
}
