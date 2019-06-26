#include <iostream>
#include <fstream>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_whitelist.h"

void bustools_whitelist(Bustools_opt &opt) {
  BUSHeader h;
  size_t nr = 0;
  size_t N = 100000;
  BUSData *p = new BUSData[N];

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

  while (true) {
    in.read((char*) p, N * sizeof(BUSData));
    size_t rc = in.gcount() / sizeof(BUSData);
    if (rc == 0) {
      break;
    }
    nr += rc;

    for (size_t i = 0; i < rc; i++) {
      if (curr_bc != p[i].barcode || bc_count == -1) {
        if (bc_count >= opt.threshold) {
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
    if (bc_count >= opt.threshold) {
      o << binaryToString(curr_bc, bclen) << "\n";
    }
  }
  /* Done reading BUS file. */

  delete[] p; p = nullptr;
  of.close();
  std::cerr << "Read in " << nr << " number of busrecords, wrote " << wl_count << " barcodes to whitelist" << std::endl;
}
