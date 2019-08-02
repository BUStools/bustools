#include <iostream>
#include <fstream>
#include <zlib.h>

#include "kseq.h"
#include "Common.hpp"
#include "BUSData.h"

#include "bustools_extract.h"

KSEQ_INIT(gzFile, gzread);

void bustools_extract(const Bustools_opt &opt) {
  BUSHeader h;
  size_t nr = 0;
  size_t N = 100000;
  BUSData *p = new BUSData[N];

  std::streambuf *buf = nullptr;
  std::ofstream of;
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
  
  std::vector<gzFile> fastq(opt.fastq.size());
  std::vector<kseq_t *> seq(opt.fastq.size());
  uint32_t iFastq = 0;
  for (int i = 0; i < opt.fastq.size(); ++i) {
    fastq[i] = gzopen(opt.fastq[i].c_str(), "r");
    seq[i] = kseq_init(fastq[i]);
    if (kseq_read(seq[i]) < 0) {
      std::cerr << "Error reading FASTQ" << std::endl;
      goto end_extract;
    }
  }

  while (true) {
    in.read((char *) p, N * sizeof(BUSData));
    size_t rc = in.gcount() / sizeof(BUSData);
    if (rc == 0) {
      break;
    }
    nr += rc;
    for (size_t i = 0; i < rc; ++i) {
      while (iFastq < p[i].flags) {
        for (const auto &s : seq) {
          if (kseq_read(s) < 0) {
            std::cerr << "Error reading FASTQ" << std::endl;
            goto end_extract;
          }
        }
        ++iFastq;
      }

      if (iFastq > p[i].flags) {
        std::cerr << "BUS file not sorted by flag" << std::endl;
        goto end_extract;
      }

      for (const auto &s : seq) {
        std::string name(s->name.s, s->name.l);
        std::string comment(s->comment.s, s->comment.l);
        std::string sequence(s->seq.s, s->seq.l);
        std::string qual(s->qual.s, s->qual.l);

        o << '@' << name << comment << '\n'
          << sequence << '\n'
          << '+' << name << comment << '\n'
          << qual << std::endl;
      }
    }
  }

  std::cout << "Read in " << nr << " BUS records" << std::endl;

end_extract:
  delete[] p;
  for (const auto &elt : seq) {
    kseq_destroy(elt);
  }
  for (const auto &elt : fastq) {
    gzclose(elt);
  }
}

