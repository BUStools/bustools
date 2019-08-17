#include <iostream>
#include <fstream>
#include <algorithm>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_project.h"

void bustools_project(Bustools_opt &opt) {
  BUSHeader h;
  std::ofstream of;

  /* Read and parse equivalence class files. */
  std::unordered_map<std::string, int32_t> txnames;
  parseTranscripts(opt.count_txp, txnames);

  std::vector<int32_t> genemap(txnames.size(), -1);
  std::unordered_map<std::string, int32_t> genenames;
  parseGenes(opt.count_genes, txnames, genemap, genenames);
  std::vector<std::string> genenamesinv(genenames.size(), "");
  for (const auto &gene : genenames) {
    genenamesinv[gene.second] = gene.first;
  }

  std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> ecmapinv;
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
  std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> geneEc2genesinv;
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
  

  /* Write gene EC matrix. */
  of.open(opt.output + ".ec");
  for (int i = 0; i < geneEc2genes.size(); ++i) {
    of << i << '\t';

    auto gene = geneEc2genes[i].begin();
    of << std::to_string(*gene);
    ++gene;
    for (gene; gene != geneEc2genes[i].end(); ++gene) {
      of << ',' << std::to_string(*gene);
    }
    of << std::endl;
#if 0
    std::string geneEc = "";
    for (const auto &gene : geneEc2genes[i]) {
      geneEc += ",";
      geneEc += std::to_string(gene);
    }
    of << i << '\t' << geneEc.substr(1) << '\n';
#endif
  }
  of.close();

  
  /* Write gene table. */
  of.open(opt.output + ".genes.txt");
  for (const auto & gene : genenamesinv) {
    of << gene << '\n';
  }
  of.close();


  /* Open input/output BUS files and write header. */

  std::streambuf *buf = nullptr;
  if (!opt.stream_out) {
    of.open(opt.output + ".bus"); 
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
 
  /* Process and output records. */
  size_t nr = 0;
  size_t nw = 0;
  size_t N = 100000;
  BUSData *p = new BUSData[N];
  BUSData currRec;
  // Gene EC --> counts for current barcode/UMI pair
  std::unordered_map<uint32_t, uint32_t> counts;
  
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
    /* Done going through BUSdata *p. */
  }
  /* Done reading BUS file. */

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
