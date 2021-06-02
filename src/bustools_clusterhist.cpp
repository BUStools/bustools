#include <iostream>
#include <fstream>
#include <cstring>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_clusterhist.h"
#include <random>
#include <map>

void bustools_clusterhist(Bustools_opt& opt) {
	BUSHeader h;
	size_t nr = 0;
	size_t N = 100000;
	uint32_t bclen = 0;
	BUSData* p = new BUSData[N];

	// read and parse the equivelence class files

	std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> ecmapinv;
	std::vector<std::vector<int32_t>> ecmap;

	std::unordered_map<std::string, int32_t> txnames;
	parseTranscripts(opt.count_txp, txnames);
	std::vector<int32_t> genemap(txnames.size(), -1);
	std::unordered_map<std::string, int32_t> genenames;
	parseGenes(opt.count_genes, txnames, genemap, genenames);
	parseECs(opt.count_ecs, h);
	ecmap = std::move(h.ecs);
	ecmapinv.reserve(ecmap.size());
	for (int32_t ec = 0; ec < ecmap.size(); ec++) {
		ecmapinv.insert({ ecmap[ec], ec });
	}
	std::vector<std::vector<int32_t>> ec2genes;
	create_ec2genes(ecmap, genemap, ec2genes);


	std::string gene_ofn = opt.output + "genes.txt";
	std::string histDir = opt.output + "cluster_hists/";


	std::vector<BUSData> v;
	v.reserve(N);
	uint64_t current_bc = 0xFFFFFFFFFFFFFFFFULL;
	//temporary data
	std::vector<int32_t> ecs;
	std::vector<int32_t> glist;
	ecs.reserve(100);
	std::vector<int32_t> u;
	u.reserve(100);

	//Read the cluster file
	std::vector<std::string> clusterNames;
	std::unordered_map<uint64_t, size_t> bcClusters;
	{
		std::ifstream ifs(opt.cluster_input_file);
		uint32_t flag = 0;
		while (!ifs.eof()) {
			std::string bc, cluster;
			ifs >> bc >> cluster;
			if (!cluster.empty()) {
				uint64_t bcBin = stringToBinary(bc, flag);
				std::size_t clust = -1;
				auto it = std::find(clusterNames.begin(), clusterNames.end(), cluster);
				if (it != clusterNames.end()) {
					clust = (it - clusterNames.begin());
				}
				else {
					clust = clusterNames.size();
					clusterNames.push_back(cluster);
				}
				bcClusters[bcBin] = clust;
			}
		}
	}
//	for (size_t i = 0; i < 5; ++i) {
//		std::cout << "rc: " << rc << "\n";
//	}

	//Allocate histograms

	//Indexed as gene*histmax + histIndex
	size_t n_genes = genenames.size();
	const uint32_t histmax = 100;//set molecules with more than histmax copies to histmax 
	typedef std::vector<double> Histograms;
	std::vector<Histograms> clusterHistograms;
	for (std::size_t i = 0; i < clusterNames.size(); ++i) {
		clusterHistograms.push_back(std::vector<double>(n_genes * histmax, 0));
	}


	//barcodes 
	auto write_barcode = [&](const std::vector<BUSData>& v) {
		if (v.empty()) {
			return;
		}
		auto it = bcClusters.find(v[0].barcode);
		if (it == bcClusters.end()) {
			return; //This is ok, could be a cell that has been filtered in quality check
		}

		auto& hist = clusterHistograms[it->second];

		size_t n = v.size();

		for (size_t i = 0; i < n; ) {
			size_t j = i + 1;
			for (; j < n; j++) {
				if (v[i].UMI != v[j].UMI) {
					break;
				}
			}

			// v[i..j-1] share the same UMI
			ecs.resize(0);
			uint32_t counts = 0;
			for (size_t k = i; k < j; k++) {
				ecs.push_back(v[k].ec);
				counts += v[k].count;
			}

			intersect_genes_of_ecs(ecs, ec2genes, glist);
			if (glist.size() == 1) {
				//Fill in histograms for prediction.
				if (glist[0] < n_genes) { //crasches with an invalid gene file otherwise
					hist[glist[0] * histmax + std::min(counts - 1, histmax - 1)] += 1.0; //histmax-1 since histograms[g][0] is the histogram value for 1 copy and so forth
				}
				else {
					std::cerr << "Mismatch between gene file and bus file, the bus file contains gene indices that is outside the gene range!\n";
				}
			}
			i = j; // increment
		}
	};

	for (const auto& infn : opt.files) {
		std::streambuf* inbuf;
		std::ifstream inf;
		if (!opt.stream_in) {
			inf.open(infn.c_str(), std::ios::binary);
			inbuf = inf.rdbuf();
		}
		else {
			inbuf = std::cin.rdbuf();
		}
		std::istream in(inbuf);

		parseHeader(in, h);
		bclen = h.bclen;

		int rc = 0;
		while (true) {
			in.read((char*)p, N * sizeof(BUSData));
			size_t rc = in.gcount() / sizeof(BUSData);
			//std::cout << "rc: " << rc << "\n";
			nr += rc;
			if (rc == 0) {
				break;
			}


			for (size_t i = 0; i < rc; i++) {
				if (p[i].barcode != current_bc) {
					// output whatever is in v
					if (!v.empty()) {
						write_barcode(v);
					}
					v.clear();
					current_bc = p[i].barcode;
				}
				v.push_back(p[i]);

			}
		}
		if (!v.empty()) {
			write_barcode(v);
		}

		if (!opt.stream_in) {
			inf.close();
		}
	}
	delete[] p; p = nullptr;



	// write genes file
	writeGenes(gene_ofn, genenames);

	//write histogram files
	for (std::size_t f = 0; f < clusterNames.size(); ++f) {

		auto& histograms = clusterHistograms[f];
		std::string hist_ofn = histDir + clusterNames[f] + ".txt";

		std::ofstream histof;
		histof.open(hist_ofn);

		for (std::size_t g = 0; g < genenames.size(); ++g) {
			//Indexed as gene*histmax + histIndex
			std::size_t offs = g * histmax;

			//first figure out the length of the histogram, don't write that to make the file smaller
			std::size_t histEnd = histmax - 1;
			for (; histEnd != 0; --histEnd) {
				if (histograms[offs + histEnd] != 0) {
					break;
				}
			}
			for (size_t c = 0; c <= histEnd; ++c) {
				if (c != 0) {
					histof << '\t';
				}
				histof << histograms[offs + c];
			}

			histof << "\n";
		}
		histof.close();
	}


	//std::cerr << "bad counts = " << bad_count <<", rescued  =" << rescued << ", compacted = " << compacted << std::endl;

	//std::cerr << "Read in " << nr << " BUS records" << std::endl;
}
