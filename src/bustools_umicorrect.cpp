#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <queue>
#include <functional>

#include "Common.hpp"
#include "BUSData.h"
#include "bustools_umicorrect.h"
#include <set>

typedef std::vector<int32_t> GeneSet;
typedef std::set<std::size_t> Neighbors;

struct Row {
	uint64_t UMI;
	uint32_t count;
	GeneSet genes;
	std::vector<BUSData> records;

	//Algorithm temp vars, added here for speed so we don't have to realloc, copy etc
	Neighbors outNeighbors; //Neighbors to this node (directed graph)
	Neighbors inNeighbors; //Nodes this node is neighbor to (other direction of the directed graph)
	bool toBeRem = false;
	bool handled = false;
};


typedef std::vector<Row> Rows;
typedef std::set<std::size_t> IndexSet;

void GetSubgraphMembers(Rows& rows, std::size_t index, IndexSet& indices) {
	indices.insert(index);
	for (auto n : rows[index].outNeighbors) {
		if (indices.count(n) == 0) { //check that this neighbor is not added before
			//call recursively to fill in
			GetSubgraphMembers(rows, n, indices);
		}
	}
	for (auto n : rows[index].inNeighbors) {
		if (indices.count(n) == 0) { //check that this neighbor is not added before
			//call recursively to fill in
			GetSubgraphMembers(rows, n, indices);
		}
	}
}

//Finds the UMIs within the possibleCopies set that are considered copies to the index molecule
//The indices set will be reduced to 
void GetCopiesOfOrigMolecule(Rows& rows, std::size_t index, GeneSet& gs, const IndexSet& possibleCopies, IndexSet& foundCopies) {
	for (auto n : rows[index].outNeighbors) {
		if ((possibleCopies.count(n) != 0) && (foundCopies.count(n) == 0)) { //only process those within possible copies, and don't process them twice
			auto newGs = intersect(gs, rows[n].genes);
			if (!newGs.empty()) { // We add neighbors that work with the gene intersection, the others are skipped
				gs = newGs;
				foundCopies.insert(n);
				//call recursively to fill in the rest of the tree, i.e. connected neighbors to the neighbors
				GetCopiesOfOrigMolecule(rows, n, gs, possibleCopies, foundCopies);
			}
		}
	}
}

inline void MergeMolecules(Row& orig, Row& copy)
{
	orig.count += copy.count;//probably unused after this point...
	//now merge the records
	for (size_t i = 0; i < copy.records.size(); ++i) {
		copy.records[i].UMI = orig.UMI;//correct the UMI
		bool handled = false;
		auto it = orig.records.begin();
		for (; it != orig.records.end(); ++it) {
			//check if the same EC already exists, then just add to the counts
			if (it->ec == copy.records[i].ec) {
				it->count += copy.records[i].count;
				handled = true;
				break;
			}
			//check if ec is smaller than the existing ec, if so we found the insertion point
			if (it->ec > copy.records[i].ec) {
				break;
			}
		}

		if (!handled) {
			orig.records.insert(it, copy.records[i]);
		}
	}
}

//We're using the Rows variable in the algorithm for speed, so don't make it const here
void ProcessBc(Rows& rows, std::ostream& outFile, const BUSHeader& h) {
	//Divide the barcode in two halves and hash on each - we're only looking for 1 Hamming distance umis.
	//One of the halves have to match if we have Hamming distance 1.
	std::unordered_multimap<uint64_t, size_t> leftMap, rightMap;
	uint64_t rightMask = 3;
	for (uint64_t i = 1; i < h.umilen / 2; ++i) {
		rightMask |= uint64_t(3) << i * 2;
	}
	uint64_t leftMask = 0xffffffffffffffff ^ rightMask;

	//build indices
	for (size_t i = 0; i < rows.size(); ++i) {
		leftMap.emplace(rows[i].UMI & leftMask, i);
		rightMap.emplace(rows[i].UMI & rightMask, i);
	}

	auto apa = leftMap.size();

	//first calculate neighbors to all
	for (std::size_t i = 0; i < rows.size(); ++i) {
		auto leftRange = leftMap.equal_range(rows[i].UMI & leftMask);
		auto rightRange = rightMap.equal_range(rows[i].UMI & rightMask);
		//The rule is: To be considered neighbors, reads need to
		//1. Belong to the same gene
		//2. Have Hamming dist 1
		//The direction is determined as:
		//i has j as out neighbour if rows[i].counts >= rows[j].counts*2 - 1
		//and vice versa. So, if both counts are 1, they will be neighbors to each other
		//This generates a directed graph

		//we have to loop twice, once for left and once for right
		for (auto it = leftRange.first; it != leftRange.second; ++it) {
			auto j = it->second;
			int d = hamming(rows[i].UMI, rows[j].UMI, h.umilen); //it will always encounter itself, with Hamming dist 0, so check for equals 1
			if (d == 1) {
				//get the gene intersection and see if it is zero
				//if not, count this as a neighbor
				GeneSet intersection = intersect(rows[i].genes, rows[j].genes);
				if (!intersection.empty()) {
					if (rows[i].count >= rows[j].count * 2 - 1) {
						rows[i].outNeighbors.insert(j);
						rows[j].inNeighbors.insert(i);
					}
					if (rows[j].count >= rows[i].count * 2 - 1) {
						rows[j].outNeighbors.insert(i);
						rows[i].inNeighbors.insert(j);
					}
				}
			}
		}
		for (auto it = rightRange.first; it != rightRange.second; ++it) {
			auto j = it->second;
			int d = hamming(rows[i].UMI, rows[j].UMI, h.umilen); //it will always encounter itself, with Hamming dist 0, so check for equals 1
			if (d == 1) {
				//get the gene intersection and see if it is zero
				//if not, count this as a neighbor
				GeneSet intersection = intersect(rows[i].genes, rows[j].genes);
				if (!intersection.empty()) {
					if (rows[i].count >= rows[j].count * 2 - 1) {
						rows[i].outNeighbors.insert(j);
						rows[j].inNeighbors.insert(i);
					}
					if (rows[j].count >= rows[i].count * 2 - 1) {
						rows[j].outNeighbors.insert(i);
						rows[i].inNeighbors.insert(j);
					}
				}
			}
		}
	}
	//now, execute the UMI correction method
	for (std::size_t i = 0; i < rows.size(); ++i) {
		//Check that this row has not already been handled as part of another subgraph
		if (!rows[i].handled) {
			//assemble the subgraph that this row belongs to
			std::set<std::size_t> indices;
			GetSubgraphMembers(rows, i, indices);

			//There is a while here, since in rare cases we will not be able to resolve the algorithm in the first round
			//due to that the gene intersection must not be empty in the whole subgraph that is collapsed
			while (indices.size() > 1) { //no need to process if the molecule is alone
				//So, find all nodes with no inNeighbours, these should be kept as original molecules
				//The rest of the copies should be merged if the gene intersection does not become empty. 
				//If no original molecules are found, the first
				//is simply selected as original molecule (only happens if all UMIs have 1 copy only)
				//keep only the one with the highest number of counts, the rest are considered copies
				//Then, just go through the original molecules in order and assign all neighbors to them.
				//This may in rare cases lead to that copies are assigned to the "wrong" molecule, but this
				//does not matter much: Example: 6-1-1-4 If we start with the UMI with 6 copies, both ones will
				//be assigned to it, although a better guess would be to assign one of them to 4.

				//find original molecules
				std::vector<std::size_t> originalMolecules;
				for (auto x : indices) {
					rows[x].handled = true;//fill in this at the same time so we don't look at this index again
					if (rows[x].inNeighbors.empty()) {
						originalMolecules.push_back(x);
					}
				}
				if (originalMolecules.empty()) {
					originalMolecules.push_back(*indices.begin());
				}

				//Then assign copies to original molecules
				for (auto om : originalMolecules) {
					auto genes = rows[om].genes;
					IndexSet foundCopies;
					indices.erase(om);//The original molecule is being handled now, so we don't need to handle it anymore
					GetCopiesOfOrigMolecule(rows, om, genes, indices, foundCopies);
					rows[om].genes = genes;
					for (auto i : foundCopies) {
						MergeMolecules(rows[om], rows[i]);
						rows[i].toBeRem = true; //discard the copy later - it's bus records have now been transferred to the other molecule
						indices.erase(i); //now handled, don't try to use it again in any of the loops
					}
				}
			}
			rows[i].handled = true;//in case of no neighbors (not really needed, right?)
		}
	}

	//now, write the bus records to file
	for (auto row : rows) {
		if (!row.toBeRem) {
			for (auto rec : row.records) {
				outFile.write((char*)&rec, sizeof(rec));
			}
		}
	}
}

void bustools_umicorrect(const Bustools_opt& opt) {
	BUSHeader h;
	size_t nr = 0;
	size_t N = 100000;
	size_t molecules = 10000;
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


	std::streambuf* buf = nullptr;
	std::ofstream of;

	if (!opt.stream_out) {
		of.open(opt.output, std::ios::binary);
		buf = of.rdbuf();
	}
	else {
		buf = std::cout.rdbuf();
	}
	std::ostream o(buf);


	std::vector<BUSData> v;
	v.reserve(N);

	Rows rows;//vector of Row
	rows.reserve(molecules);


	uint64_t current_bc = 0xFFFFFFFFFFFFFFFFULL;
	uint64_t current_umi = 0xFFFFFFFFFFFFFFFFULL;
	//temporary data
	std::vector<int32_t> ecs;
	std::vector<int32_t> glist;
	ecs.reserve(100);
	std::vector<int32_t> u;
	u.reserve(100);
	std::vector<int32_t> column_v;
	std::vector<std::pair<int32_t, double>> column_vp;
	
	column_vp.reserve(N);
	glist.reserve(100);

	int bad_count = 0;
	int compacted = 0;
	int rescued = 0;
	bool headerWritten = false;
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
		if (!headerWritten) {
			writeHeader(o, h);//just use the same header
			headerWritten = true;
		}
		bclen = h.bclen;

		Row row;

		int rc = 0;
		while (true) {
			in.read((char*)p, N * sizeof(BUSData));
			size_t rc = in.gcount() / sizeof(BUSData);
			nr += rc;
			if (rc == 0) {
				break;
			}


			for (size_t i = 0; i < rc; i++) {
				if (p[i].UMI != current_umi || p[i].barcode != current_bc) {
					if (!row.records.empty()) { //if for the very first round
						row.genes.clear();
						intersect_genes_of_ecs(ecs, ec2genes, row.genes);
						rows.push_back(row);
					}
					ecs.clear();
					ecs.push_back(p[i].ec);
					current_umi = p[i].UMI;
					row.records.clear();
					row.UMI = current_umi;
					row.count = 0;
				}
				if (p[i].barcode != current_bc) {
					// output whatever is in v
					if (!rows.empty()) {
						ProcessBc(rows, o, h);
					}
					rows.clear();
					current_bc = p[i].barcode;
				}

				ecs.push_back(p[i].ec);
				row.count += p[i].count;
				row.records.push_back(p[i]);

			}
		}
		if (!row.records.empty()) { //if for the very first round
			row.genes.clear();
			intersect_genes_of_ecs(ecs, ec2genes, row.genes);
			rows.push_back(row);
		}

		if (!rows.empty()) {
			ProcessBc(rows, o, h);
		}

		if (!opt.stream_in) {
			inf.close();
		}
	}
	delete[] p; p = nullptr;

	of.close();

	std::cerr << "Read in " << nr << " BUS records" << std::endl;
}
