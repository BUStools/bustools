#include <iostream>
#include <fstream>
#include <algorithm>
#include <queue>
#include <functional>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_merge.h"

#define TP std::pair<BUSData, int>

inline bool ncmp(const TP &a, const TP &b) {
  if (a.first.flags == b.first.flags) {
    if (a.first.barcode == b.first.barcode) {
      if (a.first.UMI == b.first.UMI) {
        if (a.first.ec == b.first.ec) {
          return a.second > b.second;
        } else {
          return a.first.ec > b.first.ec;
        }
      } else {
        return a.first.UMI > b.first.UMI;
      }
    } else {
      return a.first.barcode > b.first.barcode;
    }
  } else {
    return a.first.flags > b.first.flags;
  }
}

void bustools_merge(const Bustools_opt &opt) {
  int k = opt.files.size();
  std::vector<std::ifstream> bf(k);
  std::vector<BUSHeader> vh;
 
  /* Parse all headers. */ 
  // TODO: check for compatible headers, version numbers umi and bclen
  for (int i = 0; i < k; ++i) {
    bf[i].open((opt.files[i] + "/output.bus").c_str(), std::ios::binary);
    BUSHeader h;
    parseHeader(bf[i], h);
    
    parseECs(opt.files[i] + "/matrix.ec", h);
    vh.push_back(std::move(h));
  }

  /* Create master ec(?) */
  BUSHeader oh;
  oh.version = BUSFORMAT_VERSION;
  oh.text = "Merged files from BUStools";
  //TODO: parse the transcripts file, check that they are identical and merge.
  oh.bclen = vh[0].bclen;
  oh.umilen = vh[0].umilen;
  std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> ecmapinv;
  std::vector<std::vector<int32_t>> ectrans;        
  std::vector<int32_t> ctrans;
  
  oh.ecs = vh[0].ecs; // copy operator

  for (int32_t ec = 0; ec < oh.ecs.size(); ec++) {
    ctrans.push_back(ec);
    const auto &v = oh.ecs[ec];
    ecmapinv.insert({v, ec});
  }
  ectrans.push_back(std::move(ctrans));
  
  for (int i = 1; i < k; i++) {
    ctrans.clear();
    // merge the rest of the ecs
    int j = -1;
    for (const auto &v : vh[i].ecs) {
      j++;
      int32_t ec = -1;
      auto it = ecmapinv.find(v);
      if (it != ecmapinv.end()) {
        ec = it->second;              
      } else {
        ec = ecmapinv.size();
        oh.ecs.push_back(v); // copy
        ecmapinv.insert({v,ec});
      }
      ctrans.push_back(ec);
    }
    ectrans.push_back(ctrans);
  }

  std::vector<std::vector<int32_t>> ecmap(ecmapinv.size());
  for (const auto &ec : ecmapinv) {
    ecmap[ec.second] = ec.first;
  }
 
  /* Process data. */
  std::ofstream outf(opt.output + "/output.bus");
  writeHeader(outf, oh);

  size_t nr = 0, nw = 0;
  std::priority_queue<TP, std::vector<TP>, std::function<bool(const TP &a, const TP &b)>> pq(ncmp);
  BUSData t;
  for (int i = 0; i < k; ++i) {
    bf[i].read((char *) &t, sizeof(t));
    pq.push({t, i});
    ++nr;
  }

  BUSData curr = pq.top().first;
  curr.count = 0; // We'll count this again in the first loop
  std::unordered_set<int32_t> currec;
  while (!pq.empty()) {
    TP min = pq.top();
    pq.pop();

    BUSData &m = min.first;
    int i = min.second;
    // Do I have to check the other fields?
    if (m.flags == curr.flags && m.barcode == curr.barcode && m.UMI == curr.UMI) {
      // Same data, increase count
      curr.count += m.count;
      currec.insert(ectrans[i][m.ec]);
    } else {
      // Create new ec if necessary
      if (currec.size() == 1) {
        curr.ec = *currec.begin();
      } else {
        std::vector<int32_t> tx;
        for (const auto &ec : currec) {
          const auto &v = ecmap[ec];
          tx.insert(tx.end(), v.begin(), v.end());
        }
        std::sort(tx.begin(), tx.end());
        tx.erase(std::unique(tx.begin(), tx.end()), tx.end());

        auto it = ecmapinv.find(tx);
        if (it == ecmapinv.end()) {
          curr.ec = ecmapinv.size();
          oh.ecs.push_back(tx); // Copy
          ecmapinv.insert({tx, curr.ec});
          ecmap.push_back(tx);
        } else {
          curr.ec = it->second;
        }
      }

      if (curr.count > 0) {
        outf.write((char *) &curr, sizeof(curr));
        ++nw;
      }
      curr = m;
      currec.clear();
      currec.insert(ectrans[i][m.ec]);
    }

    // Read next
    if (bf[i].good()) {
      bf[i].read((char *) &t, sizeof(t));
      if (bf[i].gcount() > 0) {
        pq.push({t, i});
        ++nr;
      }
    }
  }

  // Write out remaining straggler
  if (curr.count > 0) {
    outf.write((char *) &curr, sizeof(curr));
    ++nw;
  }

  for (int i = 0; i < k; ++i) {
    bf[i].close();
  }
  
  /* Master ec file. */
  writeECs(opt.output + "/matrix.ec", oh);

  std::cout << "Read in " << nr << " BUS records, wrote " << nw << " BUS records" << std::endl;
}

