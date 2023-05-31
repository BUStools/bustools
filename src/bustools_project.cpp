#include <iostream>
#include <fstream>
#include <algorithm>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_project.h"

void bustools_project(Bustools_opt &opt) {
  uint32_t bclen = 0; 
  uint32_t wc_bclen = 0;
  uint32_t umilen = 0;
  BUSHeader h;
  size_t nr = 0;
  size_t N = 100000;
  BUSData* p = new BUSData[N];
  char magic[4];     
  uint32_t version = 0;
  size_t stat_map = 0;
  size_t stat_unmap = 0;

  u_map_<uint64_t, uint64_t> project_map;

  /* Load the map into project_map variable
  parse bus records and map each object (barcode, umi) with project_map
  */
  
  if (opt.type == PROJECT_BC) {
    parse_ProjectMap(opt.map, project_map); // project map is an unordered map

    std::cerr << "Found " << project_map.size() << " elements in map" << std::endl;
    std::streambuf *buf = nullptr;
    std::ofstream busf_out;

    if (!opt.stream_out) {
      busf_out.open(opt.output , std::ios::out | std::ios::binary);
      buf = busf_out.rdbuf();
    } else {
      buf = std::cout.rdbuf();
    }
    std::ostream bus_out(buf);

    bool outheader_written = false;



    nr = 0;
    BUSData bd;
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

      if (!outheader_written) {
        writeHeader(bus_out, h);
        outheader_written = true;
      }

      uint32_t bclen = h.bclen;
      uint32_t umilen = h.umilen;
      int rc = 0;
      while (true) {
        in.read((char*)p, N*sizeof(BUSData));
        size_t rc = in.gcount() / sizeof(BUSData);
        if (rc == 0) {
          break;
        }
        nr += rc;
        
        for (size_t i = 0; i < rc; i++) {
          bd = p[i];
          auto it = project_map.find(bd.barcode);
          // if the source object is in the map then replace it
          if (it != project_map.end()){
            stat_map++;
            uint64_t dest = project_map[bd.barcode];
            bd.barcode = dest;
            bus_out.write((char*) &bd, sizeof(bd));
          } else { // just write the source
            stat_unmap++;
            bus_out.write((char*) &bd, sizeof(bd));
          }

        }
      }
    }
    delete[] p; p = nullptr;
    if (!opt.stream_out) {
      busf_out.close();
    }
    std::cerr << "Processed " << nr << " BUS records" << std::endl
    << "Mapped   = " << stat_map << std::endl
    << "Unmapped = " << stat_unmap << std::endl;
  }
  if (opt.type == PROJECT_UMI) {
    parse_ProjectMap(opt.map, project_map); // project map is an unordered map

    std::cout << "Found " << project_map.size() << " elements in map" << std::endl;
    std::streambuf *buf = nullptr;
    std::ofstream busf_out;

    if (!opt.stream_out) {
      busf_out.open(opt.output , std::ios::out | std::ios::binary);
      buf = busf_out.rdbuf();
    } else {
      buf = std::cout.rdbuf();
    }
    std::ostream bus_out(buf);

    bool outheader_written = false;



    nr = 0;
    BUSData bd;
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

      if (!outheader_written) {
        writeHeader(bus_out, h);
        outheader_written = true;
      }

      uint32_t bclen = h.bclen;
      uint32_t umilen = h.umilen;
      int rc = 0;
      while (true) {
        in.read((char*)p, N*sizeof(BUSData));
        size_t rc = in.gcount() / sizeof(BUSData);
        if (rc == 0) {
          break;
        }
        nr += rc;
        
        for (size_t i = 0; i < rc; i++) {
          bd = p[i];
          auto it = project_map.find(bd.UMI);
          if (it != project_map.end()){
            stat_map++;
            uint64_t dest = project_map[bd.UMI];
            bd.UMI = dest;
            bus_out.write((char*) &bd, sizeof(bd));
          } else {
            stat_unmap++;
            bus_out.write((char*) &bd, sizeof(bd));
          }

        }
      }
    }
    delete[] p; p = nullptr;
    if (!opt.stream_out) {
      busf_out.close();
    }
    std::cerr << "Processed " << nr << " BUS records" << std::endl
    << "Mapped   = " << stat_map << std::endl
    << "Unmapped = " << stat_unmap << std::endl;
  }
  if (opt.type == PROJECT_TX) {
    std::ofstream of;
    u_map_<std::string, int32_t> txnames;
    parseTranscripts(opt.count_txp, txnames);

    std::vector<int32_t> genemap(txnames.size(), -1);
    u_map_<std::string, int32_t> genenames;
    parseGenes(opt.map, txnames, genemap, genenames);
    std::vector<std::string> genenamesinv(genenames.size(), "");
    for (const auto &gene : genenames) {
      genenamesinv[gene.second] = gene.first;
    }

    u_map_<std::vector<int32_t>, int32_t, SortedVectorHasher> ecmapinv;
    std::vector<std::vector<int32_t>> ecmap;
    parseECs(opt.count_ecs, h);
    ecmap = std::move(h.ecs);
    ecmapinv.reserve(ecmap.size());
    for (int32_t ec = 0; ec < ecmap.size(); ec++) {
      ecmapinv.insert({ecmap[ec], ec});
    }

    std::vector<std::vector<int32_t>> ec2genes;        
    create_ec2genes(ecmap, genemap, ec2genes);

    std::vector<std::vector<int32_t>> geneEc2genes = ec2genes;
    u_map_<std::vector<int32_t>, int32_t, SortedVectorHasher> geneEc2genesinv;
    std::sort(geneEc2genes.begin(), geneEc2genes.end());
    auto firstNonempty = geneEc2genes.begin();
    while (firstNonempty->size() == 0 && firstNonempty != geneEc2genes.end()) {
      ++firstNonempty;
    }
    geneEc2genes.erase(geneEc2genes.begin(), firstNonempty);
    geneEc2genes.erase(std::unique(geneEc2genes.begin(), geneEc2genes.end()), geneEc2genes.end());
    for (int32_t ec = 0; ec < geneEc2genes.size(); ++ec) {
      geneEc2genesinv.insert({geneEc2genes[ec], ec});
    }
    std::vector<int32_t> v;
    geneEc2genesinv.insert({std::move(v), -1});

    std::vector<int32_t> txEc2geneEc;
    for (const auto &txEc : ec2genes) {
      txEc2geneEc.push_back(geneEc2genesinv.at(txEc));
    }
    

    // Write gene EC matrix.
    of.open(opt.output_folder + "/matrix.ec");
    for (int i = 0; i < geneEc2genes.size(); ++i) {
      of << i << '\t';

      int first = true;
      for (const auto &gene : geneEc2genes[i]) {
        if (first) {
          first = false;
        } else {
          of << ',';
        }
        of << std::to_string(gene);
      }
      of << '\n';
    }
    of.close();

    
    // Write gene table. 
    of.open(opt.output_folder + "/genes.txt");
    for (const auto & gene : genenamesinv) {
      of << gene << '\n';
    }
    of.close();


    // Open input/output BUS files and write header.

    std::streambuf *buf = nullptr;
    if (!opt.stream_out) {
      of.open(opt.output); 
      buf = of.rdbuf();
    } else {
      buf = std::cout.rdbuf();
    }
    std::ostream o(buf);

    std::streambuf *inbuf;
    std::ifstream inf;
    if (!opt.stream_in) {
      inf.open(opt.files[0].c_str(), std::ios::binary);
      inbuf = inf.rdbuf();
    } else {
      inbuf = std::cin.rdbuf();
    }
    std::istream in(inbuf);
    parseHeader(in, h);

    h.transcripts.clear();
    for (const auto & gene : genenamesinv) {
      h.transcripts.emplace_back(gene);
    }

    h.ecs = std::move(geneEc2genes);

    writeHeader(o, h);
  
    // Process and output records.
    size_t nr = 0;
    size_t nw = 0;
    size_t N = 100000;
    BUSData *p = new BUSData[N];
    BUSData currRec;
    // Gene EC --> counts for current barcode/UMI pair
    u_map_<uint32_t, uint32_t> counts;
    
    while (true) {
      in.read((char*) p, N * sizeof(BUSData));
      size_t rc = in.gcount() / sizeof(BUSData);
      if (rc == 0) {
        break;
      }
      nr += rc;

      for (size_t i = 0; i < rc; i++) {
        if (currRec.barcode != p[i].barcode || currRec.UMI != p[i].UMI) {
          // Output BUG record
          for (const auto &rec : counts) {
            ++nw;
            currRec.ec = rec.first;
            currRec.count = rec.second;
            o.write((char *) &currRec, sizeof(BUSData));
          }
          
          currRec.barcode = p[i].barcode;
          currRec.UMI = p[i].UMI;
          counts.clear();
        }
        // Get gene EC and add entry to map
        int32_t geneEc = txEc2geneEc[p[i].ec];
        if (geneEc != -1) {
          auto it = counts.find(geneEc);
          if (it == counts.end()) {
            counts.insert({geneEc, p[i].count});
          } else {
            it->second += p[i].count;
          }
        }
      }
      // Done going through BUSdata *p.
    }
    // Done reading BUS file. 

    for (const auto &rec : counts) {
      ++nw;
      currRec.ec = rec.first;
      currRec.count = rec.second;
      o.write((char *) &currRec, sizeof(BUSData));
    }

    delete[] p; p = nullptr;
    of.close();


    std::cerr << "Read in " << nr << " BUS records, wrote " << nw << " BUG records" << std::endl;
  }
/////////////////////////////////////////////////////////////////////
/*
   if (opt.type == PROJECT_TX) { 
    
  }
  */
  
}
