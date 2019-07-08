#include <iostream>
#include <fstream>
#include <cstring>
#include <map>
#include <time.h>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_linker.h"

void bustools_linker(Bustools_opt &opt) {
  BUSHeader h;

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

  int start = opt.start, end = opt.end;
  if (start == -1) {
    start = 0;
  }
  if (end == -1) {
    end = h.bclen;
  }
  if (start >= end) {
    std::cerr << "ERROR: start or end longer than barcode length of " << std::to_string(h.bclen) << std::endl;
    exit(1);
  }
  if (start == 0 && end == h.bclen) {
    std::cerr << "ERROR: given coordinates remove entire barcode" << std::endl;
    exit(1);
  }

  uint64_t preMask = (-1) << (2 * (h.bclen - start));
  uint64_t sufMask = (uint64_t) (-1) >> (2 * (32 - h.bclen + end));
  int preShift = 2 * (end - start);

  std::streambuf *buf = nullptr;
  std::ofstream of;
  if (!opt.stream_out) {
    of.open(opt.output); 
    buf = of.rdbuf();
  } else {
    buf = std::cout.rdbuf();
  }
  std::ostream o(buf);

  h.bclen -= end - start;
  writeHeader(o, h);

  size_t nr = 0;
  size_t nw = 0;
  size_t N = 100000;
  BUSData *p = new BUSData[N];

  while (true) {
    in.read((char*) p, N * sizeof(BUSData));
    size_t rc = in.gcount() / sizeof(BUSData);
    if (rc == 0) {
      break;
    }
    nr += rc;

    for (size_t i = 0; i < rc; i++) {
      uint64_t prefix = p[i].barcode & preMask;
      prefix >>= preShift;
      uint64_t suffix = p[i].barcode & sufMask;
      p[i].barcode = prefix + suffix;
      o.write((char *) &p[i], sizeof(BUSData));
      ++nw;
    }
    /* Done going through BUSdata *p. */

  }
  /* Done reading BUS file. */
  
  delete[] p; p = nullptr;
  of.close();
  std::cerr << "Read in " << nr << " BUS records, wrote " << nw << " BUS records" << std::endl;
}

