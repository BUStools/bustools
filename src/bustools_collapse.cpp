#include <iostream>
#include <fstream>
#include <cstring>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_collapse.h"


void bustools_collapse(Bustools_opt &opt) {
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
    ecmapinv.insert({ecmap[ec], ec});
  }
  std::vector<std::vector<int32_t>> ec2genes;        
  create_ec2genes(ecmap, genemap, ec2genes);


  std::string bug_ofn = opt.output + ".bus";
  std::string gene_ofn = opt.output + ".genes.txt";

  std::streambuf* buf = nullptr;
  std::ofstream of;

  if (!opt.stream_out) {
	  of.open(bug_ofn, std::ios::binary);
	  if (of.fail()) {
		  std::cerr << "Failed to open file for writing: " << bug_ofn << std::endl;
	  }
	  buf = of.rdbuf();
  }
  else {
	  buf = std::cout.rdbuf();
  }
  std::ostream bus_out(buf);


  std::vector<BUSData> v;
  v.reserve(N);
  uint64_t current_umi = 0xFFFFFFFFFFFFFFFFULL;
  uint64_t current_bc = 0xFFFFFFFFFFFFFFFFULL;
  //temporary data
  std::vector<int32_t> ecs;
  std::vector<int32_t> glist;
  ecs.reserve(100);
  std::vector<std::pair<int32_t, double>> column_vp;
  column_vp.reserve(N);
  glist.reserve(100);

  int discarded = 0;

  auto write_umi_bug = [&](const std::vector<BUSData> &v) {
    if(v.empty()) {
      return;
    }
    column_vp.resize(0);
    
    double val = 0.0;
    size_t n = v.size();

    // v[i..j-1] share the same UMI
    ecs.resize(0);
	int32_t copies = 0;
    for (size_t i = 0; i < n; ++i) {
		ecs.push_back(v[i].ec);
		copies += v[i].count;
    }

    //There's room for something smarter here to avoid discarding UMIs with one 
	//read pointing to a non-fitting gene... I guess read errors could lead to this?
	  
	intersect_genes_of_ecs(ecs,ec2genes, glist);
	if (glist.size() == 1) {
		auto entry = v[0];
		entry.count = copies;
		entry.ec = glist[1];
		bus_out.write((char*)&entry, sizeof(entry));
	} else {
		discarded++;
	}
  };

  bool bHeaderWritten = false;

  for (const auto& infn : opt.files) { 
    std::streambuf *inbuf;
    std::ifstream inf;
    if (!opt.stream_in) {
      inf.open(infn.c_str(), std::ios::binary);
      inbuf = inf.rdbuf();
    } else {
      inbuf = std::cin.rdbuf();
    }
    std::istream in(inbuf); 

    parseHeader(in, h);
    bclen = h.bclen;
    
	if (!bHeaderWritten) {
		bHeaderWritten = true;
		// write out the initial header
		writeHeader(bus_out, h);// TODO: Think through what to write here, i.e. if the header should be a "bug" header of some kind
	}


    int rc = 0;
    while (true) {
      in.read((char*)p, N*sizeof(BUSData));
      size_t rc = in.gcount() / sizeof(BUSData);
      nr += rc;
      if (rc == 0) {
        break;
      }

      for (size_t i = 0; i < rc; i++) {
        if (p[i].UMI != current_umi || p[i].barcode != current_bc) {
          // output whatever is in v
          if (!v.empty()) {
            write_umi_bug(v);
          }
          v.clear();
          current_umi = p[i].UMI;
		  current_bc = p[i].barcode;
        }
        v.push_back(p[i]);

      }            
    }
    if (!v.empty()) {
		write_umi_bug(v);
    }

    if (!opt.stream_in) {
      inf.close();
    }
  }
  delete[] p; p = nullptr;

  //write the genes file:
  writeGenes(gene_ofn, genenames);


  std::cerr << "Read in " << nr << " BUS records, discarded " << discarded << " UMIs" << std::endl;
}
