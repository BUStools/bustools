#include <iostream>
#include <fstream>
#include <algorithm>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_whitelist.h"

#define ERROR_RATE 0.01

void bustools_whitelist(Bustools_opt &opt) {
  BUSHeader h;
  size_t nr = 0;
  size_t N = 100000;
  BUSData *p = new BUSData[N];

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
  size_t rc = 1; // Non-zero so that second while loop works when using custom threshold
  int threshold;
  int wl_count;
  
  uint32_t bc_r = 0, bc_u = 0;
  uint64_t curr_umi;
  uint64_t curr_bc;
  int bc_count = -1;

 /* Determine threshold. */ 
  if (opt.threshold) { // Custom threshold
    threshold = opt.threshold;
  } else { // Determine threshold from BUS data
    /* Get counts for all barcodes in first >=200 barcodes. */
    std::vector<wl_Record> vec;

    while (true) { 
      in.read((char*) p, N * sizeof(BUSData));
      rc = in.gcount() / sizeof(BUSData);
      if (rc == 0) {
        break;
      }
      nr += rc;

      for (size_t i = 0; i < rc; i++) {
        if (curr_bc != p[i].barcode) {
          if (bc_count != -1) {
            vec.emplace_back(curr_bc, bc_r, bc_u, bc_count);
          }
          curr_bc = p[i].barcode;
          bc_r = 1;
          bc_u = 1;
          bc_count = p[i].count;
        } else {
          ++bc_r;
          if (curr_umi != p[i].UMI) {
            if (bc_u == -1) {
              bc_u = 1;
              curr_umi = p[i].UMI;
            } else {
              ++bc_u;
            }
          }
          bc_count += p[i].count;
        }
      }
      /* Done going through BUSdata *p. */

      if (bc_count != -1) {
        vec.emplace_back(curr_bc, bc_r, bc_u, bc_count);
      }

      if (vec.size() >= 200) {
        break;
      }
    }
    /* Done retrieving first 200 barcodes. */
    // Note that the last-seen barcode may not have been fully processed by this point

    /* Sort. */
    std::sort(vec.begin(), vec.end(), [&](const wl_Record &a, const wl_Record &b) {
          if (a.count == b.count) {
            return a.barcode < b.barcode;
          } else {
            return a.count > b.count;
          }
        }
    );

    /* Determine threshold. */
    int M = 10; // Use first 10 barcodes
    int avgCount = 0, avgR = 0, avgU = 0;
    for (int i = 0; i < M; ++i) {
      avgCount += vec[i].count;
    }
    avgCount /= M;
    // [average count of top 10] * [chance of perfect barcode]
    // = [expected number of perfect barcodes]
    // And then multiply by some constant(?)
    threshold = avgCount * (1 - pow(1 - ERROR_RATE, bclen));
  
    /* Process all the records we just went through. */
    for (const auto &rec : vec) {
      if (rec.count >= threshold) {
        o << binaryToString(rec.barcode, bclen) << "\n";
        ++wl_count;
      }
    }
  }
  /* Done determining threshold. */


  /* Go through remainder of records. */
  while (rc) {
    in.read((char*) p, N * sizeof(BUSData));
    rc = in.gcount() / sizeof(BUSData);
    if (rc == 0) {
      break;
    }
    nr += rc;

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

  }
  /* Done reading BUS file. */
  
  if (bc_count >= threshold) {
    o << binaryToString(curr_bc, bclen) << "\n";
  }

  delete[] p; p = nullptr;
  of.close();
  std::cerr << "Read in " << nr << " BUS records, wrote " << wl_count << " barcodes to whitelist with threshold " << threshold << std::endl;
}
