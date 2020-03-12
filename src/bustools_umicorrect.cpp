#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <queue>
#include <functional>

#include "Common.hpp"
#include "BUSData.h"
#include "bustools_text.h"
#include <set>

typedef std::vector<int32_t> GeneSet;
typedef std::set<std::size_t> Neighbors;

struct Row {
	BUSData bd;//ec in here not used
	GeneSet genes;
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

//We're using the Rows variable in the algorithm for speed, so don't make it const here
void ProcessBc(Rows& rows, std::ostream& outFile, const BUSHeader& h) {
	//first calculate neighbors to all
	for (std::size_t i = 0; i < rows.size(); ++i) {
		for (std::size_t j = i + 1; j < rows.size(); ++j) {
			int d = hamming(rows[i].bd.UMI, rows[j].bd.UMI, h.umilen);
			//The rule is: To be considered neighbors, reads need to
			//1. Belong to the same gene
			//2. Have Hamming dist 1
			//The direction is determined as:
			//i has j as out neighbour if rows[i].counts >= rows[j].counts*2 - 1
			//and vice versa. So, if both counts are 1, they will be neighbors to each other
			//This generates a directed graph
			if (d == 1) {
				//get the gene intersection and see if it is zero
				//if not, count this as a neighbor
				GeneSet intersection = intersect(rows[i].genes, rows[j].genes);
				if (!intersection.empty()) {
					if (rows[i].bd.count >= rows[j].bd.count * 2 - 1) {
						rows[i].outNeighbors.insert(j);
						rows[j].inNeighbors.insert(i);
					}
					if (rows[j].bd.count >= rows[i].bd.count * 2 - 1) {
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
						rows[om].bd.count += rows[i].bd.count;
						rows[i].toBeRem = true;
						indices.erase(i); //now handled, don't try to use it again in any of the loops
					}
				}
			}
			rows[i].handled = true;//in case of no neighbors (not really needed, right?)
		}
	}

	//now, remove the ones marked for deletion
	for (auto it = rows.begin(); it != rows.end();) {
		if (it->toBeRem) {
			it = rows.erase(it);
		}
		else {
			++it;
		}
	}

	//Now, write the file
	for (auto row : rows) {
		write_bug_entry(outFile, row.bd, &row.genes);
	}
}

void bustools_umicorrect(const Bustools_opt& opt) {
	BUSHeader h;
	size_t nr = 0;
	size_t N = 100000;
	BUSData* p = new BUSData[N];

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

	bool bHeaderWritten = false;

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

		if (!bHeaderWritten) {
			bHeaderWritten = true;
			// write out the initial header
			writeHeader(o, h);
		}

		int rc = 0;
		size_t leftToRead = 0;
		uint64_t current_bc = 0xFFFFFFFFFFFFFFFFULL;
		Row row;
		Rows rows;

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
							ProcessBc(rows, o, h);
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
							ProcessBc(rows, o, h);
							rows.clear();
							current_bc = row.bd.barcode;
						}
						rows.push_back(row);//some unnecessary copying here, fix later
					}
				}
			}
		}
		ProcessBc(rows, o, h);
	}
	delete[] p; p = nullptr;
	if (!opt.stream_out) {
		of.close();
	}
	std::cerr << "Read in " << nr << " BUS records" << std::endl;
}
