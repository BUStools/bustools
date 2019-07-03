#include <iostream>
#include <fstream>
#include <cstring>
#include <time.h>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_inspect.h"


void bustools_inspect(Bustools_opt &opt) {
  time_t startTime = time(0);

  BUSHeader h;

  /* Load matrix.ec. */
  std::vector<std::vector<int32_t>> ecmap;
  if (opt.count_ecs.size()) {
    parseECs(opt.count_ecs, h);
    ecmap = std::move(h.ecs);
  }
  int32_t numTargets = 0;
  for (const auto &ec : ecmap) {
    if (ec.size() == 1) {
      ++numTargets;
    }
  }

  std::cerr << time(0) - startTime << ": Loaded matrix.ec" << std::endl;

  /* Load whitelist. */
  std::unordered_set<uint64_t> whitelist;
  if (opt.whitelist.size()) {
    std::ifstream wl(opt.whitelist);
    std::string inp;
    uint32_t flag;
    while (std::getline(wl, inp)) {
      whitelist.insert(stringToBinary(inp, flag));
    }
    wl.close();
  }
  
  std::cerr << time(0) - startTime << ": Loaded whitelist" << std::endl;

  /* Inspect. */
  size_t N = 100000;
  BUSData *p = new BUSData[N];

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

  /* Number of records. */
  size_t nr = 0;
  /* Number of records with multiplicity. */
  size_t reads = 0;

  /* Number of distinct barcodes. */
  uint32_t bc_count = 0;
  /* List describing number of records per barcode, with multiplicity. */
  std::vector<uint32_t> readsPerBc;
  
  /* Set of all UMIs. */
  std::unordered_set<uint64_t> umis;
  /* Number of distinct UMI/barcode pairs. */
  uint32_t umi_count = 0;
  /* List describing number of UMIs per barcode, with multiplicity. */
  std::vector<uint32_t> umisPerBc;

  /* Frequency of number of records (for Good-Toulmin).
     Multiplicity --> frequency of multiplicity. */
  std::unordered_map<uint32_t, uint32_t> freq_records;

  /* List describing number of targets per set, with multiplicity. */
  std::vector<uint32_t> targetsPerSet;
  /* Number of records with singleton target. */
  uint32_t singleton;
  /* Frequency of targets (for Good-Toulmin). */
  std::vector<uint32_t> freq_targets(numTargets, 0);

  /* Number of barcodes in agreement with whitelist. */
  uint32_t bc_wl;
  /* Number of reads in agreement with whitelist. */
  uint32_t reads_wl;

  uint64_t curr_umi, curr_bc;
  uint32_t readsPerBc_count = 0, umisPerBc_count = 0;
  bool flag = false;

  /* Process records. */
  while (true) {
    in.read((char*) p, N * sizeof(BUSData));
    size_t rc = in.gcount() / sizeof(BUSData);
    if (rc == 0) {
      break;
    }
    nr += rc;

    for (size_t i = 0; i < rc; i++) {
      if (curr_bc != p[i].barcode || !flag) {
        if (whitelist.find(curr_bc) != whitelist.end()) {
          reads_wl += readsPerBc_count;
        }

        ++bc_count;
        curr_bc = p[i].barcode;
        flag = true;
        readsPerBc.push_back(readsPerBc_count);
        readsPerBc_count = 0;

        ++umi_count;
        curr_umi = p[i].UMI; // Count distinct barcode/UMI pairs.
        umis.insert(p[i].UMI);
        umisPerBc.push_back(umisPerBc_count);
        umisPerBc_count = 1;

        if (whitelist.find(p[i].barcode) != whitelist.end()) {
          ++bc_wl;
        }
      } else if (curr_umi != p[i].UMI) {
        ++umi_count;
        curr_umi = p[i].UMI;
        umis.insert(p[i].UMI);
        ++umisPerBc_count;

      } else {
        /* Do something... */
      }

      reads += p[i].count;

      readsPerBc_count += p[i].count;

      auto ok = freq_records.insert({p[i].count, 1});
      if (!ok.second) {
        ++ok.first->second;
      }

      if (ecmap.size()) {
        for (const auto &target : ecmap[p[i].ec]) {
          ++freq_targets[target];
        }
        uint32_t targetsInSet = ecmap[p[i].ec].size();
        targetsPerSet.insert(targetsPerSet.end(), p[i].count, targetsInSet);
        if (targetsInSet == 1) {
          singleton += p[i].count;
        }
      }

    }
    /* Done going through BUSData *p. */

  }
  /* Done reading BUS file. */

  // Some stats have stragglers
  if (whitelist.find(curr_bc) != whitelist.end()) {
    reads_wl += readsPerBc_count;
  }
  readsPerBc.push_back(readsPerBc_count);
  umisPerBc.push_back(umisPerBc_count);

  delete[] p; p = nullptr;

  std::cerr << time(0) - startTime << ": Processed reads" << std::endl;

  /* Some computation. */
  size_t s;

  // Median reads per barcode
  double readsPerBcMed = 0;
  s = readsPerBc.size();
  std::sort(readsPerBc.begin(), readsPerBc.end());
  if (s > 1) { // Discard first elt = 0
    readsPerBcMed = readsPerBc[s / 2 + 1];
    if (s % 2) {
      readsPerBcMed += readsPerBc[s / 2 + 2];
      readsPerBcMed /= 2;
    }
  }
 
  // Median UMIs per barcode
  double umisPerBcMed = 0;
  s = umisPerBc.size();
  std::sort(umisPerBc.begin(), umisPerBc.end());
  if (s > 1) { // Discard first elt = 0
    umisPerBcMed = umisPerBc[s / 2 + 1];
    if (s % 2) {
      umisPerBcMed += umisPerBc[s / 2 + 2];
      umisPerBcMed /= 2;
    }
  }

  // Good-Toulmin for number of records
  int64_t gt_records = 0;
  for (const auto &elt : freq_records) {
    if (elt.first % 2) {
      gt_records += elt.second;
    } else {
      gt_records -= elt.second;
    }
  }
  
  // Median targets per set
  double targetsPerSetMed = 0;
  s = targetsPerSet.size();
  std::sort(targetsPerSet.begin(), targetsPerSet.end());
  if (s) {
    targetsPerSetMed = targetsPerSet[s / 2];
    if (s % 2) {
      targetsPerSetMed += targetsPerSet[s / 2];
      targetsPerSetMed /= 2;
    }
  }

  // Mean targets per set
  double targetsPerSetMean = 0;
  for (const auto &elt : targetsPerSet) {
    targetsPerSetMean += elt;
  }
  targetsPerSetMean /= targetsPerSet.size();

  // Good-Toulmin for number of targets
  // Also number of targets detected
  uint64_t targetsDetected = 0;
  std::unordered_map<uint32_t, uint32_t> freq_freq_targets;
  for (const auto &elt : freq_targets) {
    if (elt) {
      ++targetsDetected;
      auto ok = freq_freq_targets.insert({elt, 1});
      if (!ok.second) {
        ++ok.first->second;
      }
    }
  }
  uint64_t gt_targets = 0;
  for (const auto &elt : freq_freq_targets) {
    if (elt.first % 2) {
      gt_targets += elt.second;
    } else {
      gt_targets -= elt.second;
    }
  }
  
  std::cerr << time(0) - startTime << ": Completed calculations" << std::endl;

  /* Output info. */
  std::cout
    << "Read in " << nr << " BUS records" << std::endl
    << "Total number of reads: " << reads << std::endl
    << std::endl

    << "Number of distinct barcodes: " << std::to_string(bc_count) << std::endl
    << "Median number of reads per barcode: " << std::to_string(readsPerBcMed) << std::endl
    << "Mean number of reads per barcode: " << std::to_string((double) reads / bc_count) << std::endl
    << std::endl

    << "Number of distinct UMIs: " << std::to_string(umis.size()) << std::endl
    << "Number of distinct barcode/UMI pairs: " << std::to_string(umi_count) << std::endl
    << "Median number of UMIs per barcode: " << std::to_string(umisPerBcMed) << std::endl
    << "Mean number of UMIs per barcode: " << std::to_string((double) umi_count / bc_count) << std::endl
    << std::endl

    << "Estimated number of new records discovered if sequenced again: "
      << std::to_string(gt_records) << std::endl
    << std::endl

    << std::flush;

  if (opt.count_ecs.size()) {
    std::cout
      << "Number of distinct targets detected: " << std::to_string(targetsDetected) << std::endl
      << "Median number of targets per set: " << std::to_string(targetsPerSetMed) << std::endl
      << "Mean number of targets per set: " << std::to_string(targetsPerSetMean) << std::endl
      << std::endl

      << "Number of reads with singleton target: " << std::to_string(singleton) << std::endl
      << std::endl

      << "Estimated number of new targets discovered if sequenced again: "
        << std::to_string(gt_targets) << std::endl
      << std::endl

      << std::flush;
  }


  if (opt.whitelist.size()) {
    std::cout
      << "Number of barcodes in agreement with whitelist: " << std::to_string(bc_wl)
        << " (" << std::to_string((double) bc_wl / bc_count * 100) << "%)" << std::endl
      << "Number of reads with barcode in agreement with whitelist: " << std::to_string(reads_wl)
        << " (" << std::to_string((double) reads_wl / reads * 100) << "%)" << std::endl
      << std::endl

      << std::flush;
  }
  
  std::cerr << time(0) - startTime << ": Done" << std::endl;

}
