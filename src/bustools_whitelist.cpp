#include <iostream>
#include <fstream>
#include <algorithm>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_whitelist.h"

void bustools_whitelist(Bustools_opt &opt) {
  BUSHeader h;
  size_t nr = 0;
  size_t N = 100000;
  BUSData *p = new BUSData[N];

  int threshold;
  int wl_count;

  std::ofstream of(opt.output);
  std::ostream o(of.rdbuf());

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
  uint32_t bclen = h.bclen;

  int bc_count = -1;
  uint64_t curr_bc;
    
  in.read((char*) p, N * sizeof(BUSData));
  size_t rc = in.gcount() / sizeof(BUSData);
  nr += rc;
  
  if (opt.threshold) { // Custom threshold
    threshold = opt.threshold;
  } else { // Determine threshold from BUS data
    /* Get counts for all barcodes in first N records. */
    using Rec = std::pair<uint64_t, int>;
    std::vector<Rec> vec;
    for (size_t i = 0; i < rc; i++) {
      if (curr_bc != p[i].barcode) {
        if (bc_count != -1) {
          vec.push_back({curr_bc, bc_count});
        }

        bc_count = p[i].count;
        curr_bc = p[i].barcode;
      } else {
        bc_count += p[i].count;
      }
    }
    /* Done going through BUSdata *p. */

    if (bc_count != -1) {
      vec.push_back({curr_bc, bc_count});
    }

    /* Sort. */
    std::sort(vec.begin(), vec.end(), [&](const Rec &a, const Rec &b) {
          if (a.second == b.second) {
            return a.first < b.second;
          } else {
            return a.second > b.second;
          }
        }
    );

    /* Determine threshold. */
    int M = vec.size() / 10; // Use top 10%
    if (M < 10) {
      M = 10;
    }
    int accum = 0;
    for (int i = 0; i < M; ++i) {
      accum += vec[i].second;
    }
    // [average count of top 10%] * [chance of perfect barcode]
    // = [expected number of perfect barcodes]
    // TODO: don't hard-code the error rate, or at least make it not a magic number
    threshold = (accum / M) * pow(1 - 0.001, bclen);

    bc_count = -1; // Reset since we'll process p again
  }

  while (rc) {
    for (size_t i = 0; i < rc; i++) {
      if (curr_bc != p[i].barcode || bc_count == -1) {
        if (bc_count >= threshold) {
          o << binaryToString(curr_bc, bclen) << "\n";
          ++wl_count;
        }
        bc_count = p[i].count;
        curr_bc = p[i].barcode;
      } else {
        bc_count += p[i].count;
      }
    }
    /* Done going through BUSdata *p. */

    if (bc_count >= threshold) {
      o << binaryToString(curr_bc, bclen) << "\n";
    }
    
    in.read((char*) p, N * sizeof(BUSData));
    rc = in.gcount() / sizeof(BUSData);
    nr += rc;
  }
  /* Done reading BUS file. */

  delete[] p; p = nullptr;
  of.close();
  std::cerr << "Read in " << nr << " number of busrecords, wrote " << wl_count << " barcodes to whitelist with threshold " << threshold << std::endl;
}
