#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <queue>
#include <functional>

#include "Common.hpp"
#include "BUSData.h"
#include "bustools_text.h"
#include <algorithm>
#include <map>
#include <utility>

typedef std::vector<int32_t> GeneSet;

struct Rowq {
	BUSData bd;//ec in here not used
	GeneSet genes;
};

typedef std::vector<Rowq> Rowqs;


void bustools_predquant(const Bustools_opt& opt) {
	BUSHeader h;
	size_t nr = 0;
	size_t N = 100000;
	BUSData* p = new BUSData[N];

	std::streambuf* buf = nullptr;
	std::ofstream of;

	std::string mtx_ofn = opt.output + ".mtx";
	std::string barcodes_ofn = opt.output + ".barcodes.txt";
	std::string gene_ofn = opt.output + ".genes.txt";
	of.open(mtx_ofn, std::ios::binary); //do not change this to text mode, will mess things up for Windows since \n becomes \r\n in text mode, taking extra bytes
	buf = of.rdbuf();

	// write out the initial header
	of << "%%MatrixMarket matrix coordinate real general\n%\n";
	// number of genes
	auto mat_header_pos = of.tellp();
	std::string dummy_header(66, '\n');
	for (int i = 0; i < 33; i++) {
		dummy_header[2 * i] = '%';
	}
	of.write(dummy_header.c_str(), dummy_header.size());

	size_t n_cols = 0;
	size_t n_rows = 0;
	size_t n_entries = 0;
	size_t n_totcounts = 0;

	std::vector<uint64_t> barcodes;

	//get genes:
	std::vector<std::string> genenames;
	parseBugGenes(opt.predquant_genefile, genenames);

	n_cols = genenames.size();

	std::vector<std::pair<int32_t, double>> column_vp;
	column_vp.reserve(N);

	//Allocate histograms
	const uint32_t histmax = 100;//set molecules with more than histmax copies to histmax 
	auto histograms = new double[genenames.size()][histmax];
	memset(histograms, 0, sizeof(double) * genenames.size() * histmax);

	auto ProcessBc = [&](Rowqs& rows) {
		if (rows.empty()) {
			return;
		}
		column_vp.resize(0);
		n_rows += 1;

		barcodes.push_back(rows[0].bd.barcode);
		double val = 0.0;
		size_t n = rows.size();

		for (size_t i = 0; i < n; ++i) {
			size_t gn = rows[i].genes.size();
			for (auto x : rows[i].genes) {
				column_vp.push_back({ x, 1.0 / gn });

				//Fill in histograms for prediction later.
				//Also fill in multimapped reads for now
				if (x < n_cols) { //crasches with an invalid gene file otherwise
					histograms[x][std::min(rows[i].bd.count-1, histmax - 1)] += 1.0 / gn; //histmax-1 since histograms[g][0] is the histogram value for 1 copy and so forth
				} else {
					std::cerr << "Mismatch between gene file and bus file, the bus file contains gene indices that is outside the gene range!\n";
				}
			}

			++n_totcounts;
		}

		std::sort(column_vp.begin(), column_vp.end());
		size_t m = column_vp.size();
		std::unordered_map<int32_t, double> col_map(m);
		std::vector<int32_t> cols;

		for (size_t i = 0; i < m; ) {
			size_t j = i + 1;
			double val = column_vp[i].second;
			for (; j < m; j++) {
				if (column_vp[i].first != column_vp[j].first) {
					break;
				}
				val += column_vp[j].second;
			}
			col_map.insert({ column_vp[i].first,val });
			cols.push_back(column_vp[i].first);

			n_entries++;

			i = j; // increment
		}

		for (const auto& x : cols) {
			double val = 0;
			auto it = col_map.find(x);
			if (it != col_map.end()) {
				val = it->second;
			}
			of << n_rows << " " << (x + 1) << " " << val << "\n";
		}
	};


	for (const auto& infn : opt.files) {
		std::streambuf* inbuf;
		std::ifstream inf;
		if (!opt.stream_in) {
			inf.open(infn.c_str(), std::ios::binary);
			inbuf = inf.rdbuf();
			if (!inf) {
				std::cerr << "Failed to read file: " << infn << std::endl;
			}
		}
		else {
			inbuf = std::cin.rdbuf();
		}
		std::istream in(inbuf);

		parseHeader(in, h);
		
		int rc = 0;
		size_t leftToRead = 0;
		uint64_t current_bc = 0xFFFFFFFFFFFFFFFFULL;
		Rowq row;
		Rowqs rows;

		while (true) {
			in.read((char*)p, N * sizeof(BUSData));
			size_t rc = in.gcount() / sizeof(BUSData);
			if (rc == 0) {
				break;
			}
			nr += rc;
			for (size_t i = 0; i < rc; i++) {
				if (leftToRead > 0) {
					add_ecs_from_block(p[i], row.genes, leftToRead);
					if (leftToRead == 0) {
						if (row.bd.barcode != current_bc) {
							ProcessBc(rows);
							rows.clear();
							current_bc = row.bd.barcode;
						}
						rows.push_back(row);//some unnecessary copying here, fix later
					}
				} else {
					if (p[i].ec < 0) {
						//start a "multiblock" read
						row.genes.clear();
						row.bd = p[i];
						leftToRead = -p[i].ec;
					} else {
						row.genes.clear();
						row.genes.push_back(p[i].ec);
						row.bd = p[i];

						if (row.bd.barcode != current_bc) {
							ProcessBc(rows);
							rows.clear();
							current_bc = row.bd.barcode;
						}
						rows.push_back(row);//some unnecessary copying here, fix later
					}
				}
			}
		}
		ProcessBc(rows);
		
	}
	delete[] p; p = nullptr;
	if (!opt.stream_out) {
		of.close();
	}

	std::stringstream ss;
	ss << n_rows << " " << n_cols << " " << n_entries << "\n";
	std::string header = ss.str();
	size_t hlen = header.size();
	assert(hlen < 66);
	of.open(mtx_ofn, std::ios::binary | std::ios::in | std::ios::out);
	of.seekp(mat_header_pos);
	of.write("%", 1);
	of.write(std::string(66 - hlen - 2, ' ').c_str(), 66 - hlen - 2);
	of.write("\n", 1);
	of.write(header.c_str(), hlen);
	of.close();

	// write barcode file
	std::ofstream bcof;
	bcof.open(barcodes_ofn);
	for (const auto& x : barcodes) {
		bcof << binaryToString(x, h.bclen) << "\n";
	}
	bcof.close();

	std::cerr << "Read in " << nr << " BUS records" << std::endl;

	//now, create a scaled counts matrix taking prediction into account
	//We do Good-Toulmin for now:
	auto predict_double_good_toulmin = [&](double hist[]) {
		double predval = 0;
		for (size_t i = 0; i < histmax; i += 2) {
			predval += hist[i];
		}
		for (size_t i = 1; i < histmax; i += 2) {
			predval -= hist[i];
		}

		return predval;
	};

	auto predict = predict_double_good_toulmin; //switch to preseq based on flags later
	//if (opt.predquant_algorithm == PREDQUANT_ALG::PRESEQ;) {
	//	predict = predict_preseq, also add something about how long to predict
	//}


	//Loop through each gene and predict
	auto geneTotCounts = new double[genenames.size()];
	auto genePredCounts = new double[genenames.size()];
	auto geneFracOnes = new double[genenames.size()];
	memset(geneTotCounts, 0, sizeof(double)* genenames.size());
	memset(genePredCounts, 0, sizeof(double)* genenames.size());

	std::multimap<double, size_t> aboveBucketLevelIndices;//will automatically be sorted on bucketlimit

	//get current counts and fraction of ones
	for (size_t g = 0; g < genenames.size(); ++g) {
		double sum = 0;
		for (size_t i = 0; i < histmax; ++i) {
			sum += histograms[g][i];
		}
		geneTotCounts[g] = sum;
		geneFracOnes[g] = histograms[g][0] / sum;
		if (sum >= opt.predquant_include_bucket_limit) {
			aboveBucketLevelIndices.insert(std::make_pair(geneFracOnes[g], g));
		}
	}

	size_t numBuckets = std::min((size_t) opt.predquant_num_buckets, aboveBucketLevelIndices.size());
	auto bucketFracOnes = new double[numBuckets];
	auto bucketHistograms = new double[numBuckets][histmax];
	memset(bucketHistograms, 0, sizeof(double)*numBuckets*histmax);

	//create buckets
	size_t index = 0;
	size_t lastBucketIndex = 0;
	size_t numInBucket = 0;
	size_t totBucketGenes = aboveBucketLevelIndices.size();
	if (totBucketGenes == 0) {
		std::cerr << "Found zero genes above bucket level - need to lower the bucket inclusion threshold!\n";
		//this will likely leak memory...
		return;
	}
	double scale = double(numBuckets)/double(totBucketGenes);
	for (auto pair : aboveBucketLevelIndices) {
		//calculate bucket index
		size_t bucketIndex = size_t(double(index) * scale);
		//add histogram to bucket
		for (size_t i = 0; i < histmax; ++i) {
			bucketHistograms[bucketIndex][i] += histograms[pair.second][i];
		}
		++index;
	}
	
	//now predict the buckets, calc frac of ones and create scale factors for them
	auto bucketScaling = new double[numBuckets];
	for (size_t b = 0; b < numBuckets; ++b) {
		double sum = 0.0;
		//this loop could be optimized away if needed
		for (size_t i = 0; i < histmax; ++i) {
			sum += bucketHistograms[b][i];
		}
		bucketScaling[b] = (sum + std::max(0.0, predict(bucketHistograms[b])))/sum;
		bucketFracOnes[b] = bucketHistograms[b][0] / sum;
	}

	//we're now ready to predict all genes - the lowly expressed ones by mapping them to the closest bucket
	for (size_t g = 0; g < genenames.size(); ++g) {
		if (geneTotCounts[g] >= opt.predquant_use_bucket_limit) {
			//predict directly
			genePredCounts[g] = geneTotCounts[g] + std::max(0.0, predict(histograms[g]));
		} else {
			//map to bucket
			//loop through all buckets and see which fraction of ones that match the best
			auto fo = geneFracOnes[g];
			for (size_t i = 0; i < numBuckets; ++i) {
				//if we are at the last item in the loop, or if this item is closer than the next, then we have found our best match
				if (i == numBuckets - 1 || (fabs(fo - bucketFracOnes[i]) <= fabs(fo - bucketFracOnes[i + 1]))) {
					//scale the original counts with the scaling from the closest bucket
					genePredCounts[g] = geneTotCounts[g] * bucketScaling[i];
					break;
				}
			}
		}
	}



	//now, calculate scaling factors for each gene to apply on the counts matrix
	auto scaleFactors = new double[genenames.size()];
	//Loop through all genes, sum up the total counts and total predicted counts, 
	//and scale the predicted counts by total counts/predicted counts. We then have 
	//the scaled total counts for each gene. The scale factors are then for each gene
	//scaled predicted counts/total counts 
	//first calc sums
	double sumPredCounts = 0;
	double sumTotCounts = 0;
	for (size_t g = 0; g < genenames.size(); ++g) {
		sumPredCounts += genePredCounts[g];
		sumTotCounts += geneTotCounts[g];
	}
	//then apply on genes
	double gsf = sumTotCounts / sumPredCounts;//global scale factor
	for (size_t g = 0; g < genenames.size(); ++g) {
		scaleFactors[g] = gsf*genePredCounts[g]/ geneTotCounts[g];
	}

	//clean up, will only need the scale factors
	delete[] histograms;
	delete[] geneTotCounts;
	delete[] genePredCounts;
	delete[] geneFracOnes;
	delete[] bucketFracOnes;
	delete[] bucketHistograms;
	delete[] bucketScaling;

	//The strategy here is to now load the matrix again, it is a much smaller file than the bug file
	//and require less processing
	//Now, load the matrix file and write it to a new file name with corrected counts
	std::string mtx_outpred = opt.output + "_pred.mtx";
	//do not change these to text mode, will mess things up for Windows since \n becomes \r\n in text mode, taking extra bytes
	std::ifstream inf(mtx_ofn, std::ios::binary);//the last out file now as in file
	of.open(mtx_outpred, std::ios::binary);

	//first skip header in in file and write it in out file
	char dummybuf[1000];
	inf.read(dummybuf, mat_header_pos);
	inf.read(dummybuf, 66);
	of << "%%MatrixMarket matrix coordinate real general\n%\n";
	of.write("%", 1);
	of.write(std::string(66 - hlen - 2, ' ').c_str(), 66 - hlen - 2);
	of.write("\n", 1);
	of.write(header.c_str(), hlen);
	
	int bc = 0, gene = 0;
	double count = 0;
	while (inf >> bc >> gene >> count) {
		//don't do any boundary checks, this should be safe since we just wrote the file ourselves
		of << bc << " " << gene << " " << count * scaleFactors[gene - 1] << "\n";
	}
	
	of.close();

	
	delete[] scaleFactors;
}
