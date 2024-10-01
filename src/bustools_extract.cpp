#include <iostream>
#include <fstream>
#include <zlib.h>

#include "kseq.h"
#include "Common.hpp"
#include "BUSData.h"

#include "bustools_extract.h"

KSEQ_INIT(gzFile, gzread);

inline bool open_fastqs(
    std::vector<gzFile> &outFastq,
    std::vector<gzFile> &inFastq,
    std::vector<kseq_t *> &seq,
    const Bustools_opt &opt, size_t &iFastq) {
  
  for (int i = 0; i < opt.nFastqs; ++i) {
    gzclose(outFastq[i]);
    outFastq[i] = gzopen(std::string(opt.output + "/" + std::to_string(iFastq + 1) + ".fastq.gz").c_str(), "w");
    gzclose(inFastq[i]);
    inFastq[i] = gzopen(opt.fastq[iFastq].c_str(), "r");
  
    if (seq[i]) {
      kseq_destroy(seq[i]);
    }
    seq[i] = kseq_init(inFastq[i]);
    if (kseq_read(seq[i]) < 0) {
      return false;
    }
    
    ++iFastq;
  }
  return true;
}

void bustools_extract(const Bustools_opt &opt) {
  BUSHeader h;
  size_t nr = 0;
  size_t N = 100000;
  BUSData *p = new BUSData[N];
  char *buf = new char[N];
  buf[0] = '@';

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
  
  std::vector<gzFile> outFastq(opt.nFastqs);
  std::vector<gzFile> inFastq(opt.nFastqs);
  std::vector<kseq_t *> seq(opt.nFastqs, nullptr);
  uint32_t iRead = 0;
  size_t iFastq = 0;
  uint32_t lastFlag = 0;


  auto write_seq_to_file = [&opt, &buf, &outFastq] (std::vector<kseq_t *> &seq) {
    for (int i = 0; i < opt.nFastqs; ++i) {
        int bufLen = 1; // Already have @ character in buffer
        
        memcpy(buf + bufLen, seq[i]->name.s, seq[i]->name.l);
        bufLen += seq[i]->name.l;

        // Only add space and comment if the comment is not empty
        if (seq[i]->comment.l > 0) {
            // Add space between name and comment
            buf[bufLen++] = ' ';

            memcpy(buf + bufLen, seq[i]->comment.s, seq[i]->comment.l);
            bufLen += seq[i]->comment.l;
        }

        buf[bufLen++] = '\n';

        memcpy(buf + bufLen, seq[i]->seq.s, seq[i]->seq.l);
        bufLen += seq[i]->seq.l;
        
        buf[bufLen++] = '\n';
        buf[bufLen++] = '+';
        buf[bufLen++] = '\n';

        memcpy(buf + bufLen, seq[i]->qual.s, seq[i]->qual.l);
        bufLen += seq[i]->qual.l;
        
        buf[bufLen++] = '\n';

        if (gzwrite(outFastq[i], buf, bufLen) != bufLen) {          
          return false;
        }
      }
      return true;
  };

  bool tail = false;
  bool finished = false;
  size_t iFlag = 0;
  size_t rc = 0;
  
  if (!open_fastqs(outFastq, inFastq, seq, opt, iFastq)) {
    std::cerr << "Error reading FASTQ " << opt.fastq[iFastq] << std::endl;
    goto end_extract;
  }

  // fill in the first N BUS records
  in.read((char *) p, N * sizeof(BUSData));
  rc = in.gcount() / sizeof(BUSData);  
  nr += rc;
  tail = rc==0;

  while (true) {
    // read in the sequence    
    if (opt.extract_include && iRead == p[iFlag].flags) {
      if (!write_seq_to_file(seq)) {
        std::cerr << "Error writing to FASTQ" << std::endl;
        goto end_extract;
      }
    }
  

    if (opt.extract_exclude && (iRead < p[iFlag].flags || tail)) {
      if (!write_seq_to_file(seq)) {
        std::cerr << "Error writing to FASTQ" << std::endl;
        goto end_extract;
      }
    }

    if (!tail && iRead == p[iFlag].flags) {
      // read the next flag from the next bus record
      iFlag++;
      
      if (iFlag == rc)  {
        // read the next batch of bus
        in.read((char *) p, N * sizeof(BUSData));
        rc = in.gcount() / sizeof(BUSData);
        iFlag = 0;
        tail = rc==0;
      } 
    }
    // fill the next read
    for (const auto &s : seq) {
      int err_kseq_read = kseq_read(s);
      if (err_kseq_read == -1) { // Reached EOF
        // check if we are done with all files
        if (iFastq == opt.fastq.size()) { // Done with all files
          finished = true;
        } else {
          if (!open_fastqs(outFastq, inFastq, seq, opt, iFastq)) {
            std::cerr << "Error: cannot read FASTQ " << opt.fastq[iFastq] << std::endl;
            goto end_extract;
          }
        }
      } else if (err_kseq_read == -2) {
        std::cerr << "Error: truncated FASTQ" << std::endl;
        goto end_extract;
      }
    } 
    ++iRead;

    if (finished) {
      if (iFlag < rc) {
        std::cerr << "Warning: number of reads in FASTQs was less than number of reads in BUS file" << std::endl;
        goto end_extract;
      }
      break;
    }

    
  }

  

  std::cerr << "Read in " << nr << " BUS records" << std::endl;

end_extract:
  delete[] p;
  delete[] buf;
  for (auto &elt : outFastq) {
    gzclose(elt);
  }
  for (auto &elt : inFastq) {
    gzclose(elt);
  }
  for (auto &elt : seq) {
    if (elt) {
      kseq_destroy(elt);
    }
  }
}

