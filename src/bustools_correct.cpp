#include <iostream>
#include <fstream>
#include <algorithm>

#include <unordered_map>
#include <vector>


#include "Common.hpp"
#include "BUSData.h"

#include "bustools_correct.h"

#include "roaring.hh"


int search_for_mismatch(const Roaring& r, const size_t bc, const uint64_t b, uint64_t &c) {
  int counts = 0;
  if (r.isEmpty()) {
    return 0;
  } else {
    size_t sh = bc-1;          

    for (size_t i = 0; i < bc; ++i) {
      for (uint64_t d = 1; d <= 3; d++) {
        uint64_t y = b ^ (d << (2*sh));
        if (r.contains(y)) {
          if (counts == 0) {
            c = y;
          }
          counts++;
        }
      }                
      sh--;
    }
  }
  return counts;
}


void bustools_correct(Bustools_opt &opt) {
  uint32_t bclen = 0; 
  uint32_t wc_bclen = 0;
  uint32_t umilen = 0;
  BUSHeader h;
  size_t nr = 0;
  size_t N = 100000;
  BUSData* p = new BUSData[N];
  char magic[4];      
  uint32_t version = 0;
  size_t stat_white = 0;
  size_t stat_corr = 0;
  size_t stat_uncorr = 0;

  std::ifstream wf(opt.whitelist, std::ios::in);
  std::string line;
  line.reserve(100);
  std::unordered_set<uint64_t> wbc;
  wbc.reserve(100000);
  uint32_t f = 0;
  while(std::getline(wf, line)) {
    if (wc_bclen == 0) {
      wc_bclen = line.size();
    }
    uint64_t bc = stringToBinary(line, f);
    wbc.insert(bc);          
  }
  wf.close();

  std::cerr << "Found " << wbc.size() << " barcodes in the whitelist" << std::endl;

  // split barcode into upper and lower half
  size_t bc2 = (wc_bclen+1)/2; 

  std::vector<std::pair<Roaring,Roaring>> correct(1ULL<<(2*bc2)); // 4^(bc/2) possible barcodes

  uint64_t mask_size = (1ULL << (2*bc2));
  uint64_t lower_mask = (1ULL<<(2*bc2))-1;
  uint64_t upper_mask = (1ULL<<(2*(wc_bclen-bc2)))-1;
  for (uint64_t b : wbc) {
    uint64_t lb = b & lower_mask;
    uint64_t ub = (b>>(2*bc2)) & upper_mask;

    correct[ub].second.add(lb);
    correct[lb].first.add(ub);
  }

  std::streambuf *buf = nullptr;
  std::ofstream busf_out;
  
  if (!opt.stream_out) {
    busf_out.open(opt.output , std::ios::out | std::ios::binary);
    buf = busf_out.rdbuf();
  } else {
    buf = std::cout.rdbuf();
  }
  std::ostream bus_out(buf);

  bool outheader_written = false;
  
  nr = 0;
  BUSData bd;
  for (const auto& infn : opt.files) { 
    std::streambuf *inbuf;
    std::ifstream inf;
    if (!opt.stream_in) {
      inf.open(infn.c_str(), std::ios::binary);
      inbuf = inf.rdbuf();
    } else {
      inbuf = std::cin.rdbuf();
    }
    std::istream in(inbuf);          
    parseHeader(in, h);

    if (!outheader_written) {
      writeHeader(bus_out, h);
      outheader_written = true;
    }

    if (bclen == 0) {
      bclen = h.bclen;

      if (bclen != wc_bclen) { 
        std::cerr << "Error: barcode length and whitelist length differ, barcodes = " << bclen << ", whitelist = " << wc_bclen << std::endl
                  << "       check that your whitelist matches the technology used" << std::endl;

        exit(1);
      }
    }
    if (umilen == 0) {
      umilen = h.umilen;
    }

    int rc = 0;
    while (true) {
      in.read((char*)p, N*sizeof(BUSData));
      size_t rc = in.gcount() / sizeof(BUSData);
      if (rc == 0) {
        break;
      }
      nr +=rc;

      for (size_t i = 0; i < rc; i++) {
        bd = p[i];
        auto it = wbc.find(bd.barcode);
        if (it != wbc.end()) { 
          stat_white++;
          bd.count = 1;
          bus_out.write((char*) &bd, sizeof(bd));
        } else {  
          uint64_t b = bd.barcode;
          uint64_t lb = b & lower_mask;
          uint64_t ub = (b>>(2*bc2)) & upper_mask;   
          uint64_t lbc=0,ubc=0;
          int correct_lower = search_for_mismatch(correct[ub].second,bc2,lb,lbc);
          int correct_upper = search_for_mismatch(correct[lb].first,wc_bclen - bc2,ub,ubc);
          int nc = correct_lower + correct_upper;
          if (nc != 1) {
            stat_uncorr++;
          } else if (nc==1) {
            if (correct_lower == 1) {
              uint64_t b_corrected = (ub << (2*bc2)) | lbc;          
              bd.barcode = b_corrected;
              stat_corr++;
              
            } else if (correct_upper == 1) {
              uint64_t b_corrected = (ubc << (2*bc2)) | lb; 
              bd.barcode = b_corrected;
              stat_corr++;
            }       
            bd.count = 1;
            bus_out.write((char*) &bd, sizeof(bd));     
          }
        }
      }
    }
  }

  std::cerr << "Processed " << nr << " bus records" << std::endl
  << "In whitelist = " << stat_white << std::endl
  << "Corrected = " << stat_corr << std::endl
  << "Uncorrected = " << stat_uncorr << std::endl;


  if (!opt.stream_out) {
    busf_out.close();
  }

  delete[] p; p = nullptr;
}