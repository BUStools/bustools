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
  char *buf = new char[N];
  buf[0] = '@';

  size_t K = opt.fastq.size();

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
  
 std::vector<gzFile> of(K);
  std::vector<gzFile> fastq(K);
  std::vector<kseq_t *> seq(K, nullptr);
  uint32_t iFastq = 0;
  for (int i = 0; i < K; ++i) {
    of[i] = gzopen(std::string(opt.output + "/" + std::to_string(i + 1) + ".fastq.gz").c_str(), "w");
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

      for (int i = 0; i < K; ++i) {
        int bufLen = 1; // Already have @ character in buffer
        
        memcpy(buf + bufLen, seq[i]->name.s, seq[i]->name.l);
        bufLen += seq[i]->name.l;
        
        memcpy(buf + bufLen, seq[i]->comment.s, seq[i]->comment.l);
        bufLen += seq[i]->comment.l;
        
        buf[bufLen++] = '\n';

        memcpy(buf + bufLen, seq[i]->seq.s, seq[i]->seq.l);
        bufLen += seq[i]->seq.l;
        
        buf[bufLen++] = '\n';
        buf[bufLen++] = '+';

        memcpy(buf + bufLen, seq[i]->name.s, seq[i]->name.l);
        bufLen += seq[i]->name.l;
        
        memcpy(buf + bufLen, seq[i]->comment.s, seq[i]->comment.l);
        bufLen += seq[i]->comment.l;
        
        buf[bufLen++] = '\n';

        memcpy(buf + bufLen, seq[i]->qual.s, seq[i]->qual.l);
        bufLen += seq[i]->qual.l;
        
        buf[bufLen++] = '\n';

        if (gzwrite(of[i], buf, bufLen) != bufLen) {
          std::cerr << "Error writing to FASTQ" << std::endl;
          goto end_extract;
        }
      }
    }
  }

  std::cout << "Read in " << nr << " BUS records" << std::endl;

end_extract:
  delete[] p;
  delete[] buf;
  for (auto &elt : of) {
    gzclose(elt);
  }
  for (auto &elt : fastq) {
    gzclose(elt);
  }
  for (auto &elt : seq) {
    if (elt) {
      kseq_destroy(elt);
    }
  }
}

