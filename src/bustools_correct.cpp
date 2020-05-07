#include <iostream>
#include <fstream>
#include <algorithm>

#include <unordered_map>
#include <vector>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_correct.h"

#include "roaring.hh"

int search_for_mismatch(const Roaring &r, const size_t bc, const uint64_t b, uint64_t &c)
{
  int counts = 0;
  if (r.isEmpty())
  {
    return 0;
  }
  else
  {
    size_t sh = bc - 1;

    for (size_t i = 0; i < bc; ++i)
    {
      for (uint64_t d = 1; d <= 3; d++)
      {
        uint64_t y = b ^ (d << (2 * sh));
        if (r.contains(y))
        {
          if (counts == 0)
          {
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

bool check_distance(std::vector<std::pair<Roaring, Roaring>> &correct, uint64_t b, const size_t bc_len)
{
  // TODO: create a function that takes in a correct matrix and a partial barcode
  // and tells you true or false if it is within one hamming distance of the half-whitelist
  // Question: can I correct the first half? Can I correct the second half? Is each half
  // in the "half"-whitelist? Is the concatenated "full" barcode in the whitelist?
  return false;
}

void bustools_split_correct(Bustools_opt &opt)
{
  uint32_t bclen = 0;
  uint32_t wc_bclen = 0;
  uint32_t umilen = 0;
  BUSHeader h;
  size_t nr = 0;
  size_t N = 100000;
  BUSData *p = new BUSData[N];
  char magic[4];
  uint32_t version = 0;
  size_t stat_white = 0;
  size_t stat_corr = 0;
  size_t stat_corr_2 = 0;
  size_t stat_uncorr = 0;
  uint64_t old_barcode;

  bool dump_bool = opt.dump_bool;

  std::ofstream of;
  if (dump_bool)
  {
    of.open(opt.dump);
  }

  std::ifstream wf(opt.whitelist, std::ios::in);
  std::string line;
  line.reserve(100);

  uint32_t f = 0;

  std::unordered_set<uint64_t> wbc;
  wbc.reserve(100000);

  std::unordered_set<uint64_t> wbc_34;
  wbc_34.reserve(100000);

  std::unordered_set<uint64_t> wbc_12;
  wbc_12.reserve(100000);

  // split barcode into upper and lower half
  // If even
  // AGCC AAGA GAGG GCAT
  // UU   LU    UL  LL
  // AGCC AAGA GAGG GCAT
  // 4    3    2    1
  // if odd
  // AGCC AAGA GAGG GCA

  // note integer division discards the decimal
  // let wc_bclen = 14, 15, 16,

  //
  size_t len_12; // = wc_bclen/2 ; // = 7,7, 8
  size_t len_34; // = wc_bclen-len_12; // = 7, 8, 8

  size_t len_1; // = len_12/2; // = 3, 3, 4
  size_t len_2; // = len_12 - len_1; // =4, 4, 4
  size_t len_3; // = len_34/2; //3, 4,4
  size_t len_4; // = len_34 - len_3; //4, 4,4

  // uint64_t mask_size = (1ULL << (2 * bc2));
  uint64_t mask_12; // = (1ULL << (2 * len_12)) - 1;
  uint64_t mask_34; // = (1ULL << (2 * (len_34))) - 1;

  uint64_t mask_1 = (1ULL << (2 * len_1)) - 1;
  uint64_t mask_2 = (1ULL << (2 * len_2)) - 1;
  uint64_t mask_3 = (1ULL << (2 * len_3)) - 1;
  uint64_t mask_4 = (1ULL << (2 * len_4)) - 1;

  while (std::getline(wf, line))
  {
    if (wc_bclen == 0)
    {
      wc_bclen = line.size();

      len_12 = wc_bclen / 2;      // = 7,7, 8
      len_34 = wc_bclen - len_12; // = 7, 8, 8

      // uint64_t mask_size = (1ULL << (2 * bc2));
      mask_12 = (1ULL << (2 * len_12)) - 1;
      mask_34 = (1ULL << (2 * (len_34))) - 1;
    }
    uint64_t bc = stringToBinary(line, f);
    //wbc.insert(bc);

    uint64_t bc_12 = bc & mask_12;
    uint64_t bc_34 = (bc >> (2 * len_12)) & mask_34;

    wbc_12.insert(bc_12);
    wbc_34.insert(bc_34);

    //std::cout << line<< "\t"<<binaryToString(bc_34, len_34) << "\t" << binaryToString(bc_12, len_12) << "\n";
  }
  wf.close();

  std::cerr << "Found " << wbc_12.size() << "," << wbc_34.size() << " barcodes in the half whitelists" << std::endl;

  len_1 = len_12 / 2;     // = 3, 3, 4
  len_2 = len_12 - len_1; // =4, 4, 4
  len_3 = len_34 / 2;     //3, 4,4
  len_4 = len_34 - len_3; //4, 4,4

  mask_1 = (1ULL << (2 * len_1)) - 1;
  mask_2 = (1ULL << (2 * len_2)) - 1;
  mask_3 = (1ULL << (2 * len_3)) - 1;
  mask_4 = (1ULL << (2 * len_4)) - 1;

  // std::vector<std::pair<Roaring, Roaring>> correct(1ULL << (2 * bc2)); // 4^(bc/2) possible barcodes

  std::vector<std::pair<Roaring, Roaring>> correct_12(1ULL << (2 * len_12)); // 4^(bc/4) possible barcodes
  std::vector<std::pair<Roaring, Roaring>> correct_34(1ULL << (2 * len_34)); // 4^(bc/4) possible barcodes

  // std::vector<uint64_t> masks = {};

  // barcode is split into 4 peices of equal size, check each half in the following way:
  // AGCC AAGA GAGG GCAT
  // UU   LU    UL  LL
  // U: upper, L: lower
  // create all hamming 1 distance variants of the LL half and check to see if
  // any one of those variants (concatenated with the UL) is in the "half"-whitelist
  // Do the same for the UU and LU half.

  for (uint64_t b : wbc_12)
  {

    uint64_t bc1 = b & mask_1;
    uint64_t bc2 = (b >> (2 * len_1)) & mask_2;

    // std::cout << binaryToString(b, len_12) << "\t" << binaryToString(bc2, len_2) << "\t" << binaryToString(bc1 , len_1) << "\n";

    correct_12[bc2].second.add(bc1);
    correct_12[bc1].first.add(bc2);
  }

  for (uint64_t b : wbc_34)
  {
    uint64_t bc3 = b & mask_3;
    uint64_t bc4 = (b >> (2 * len_3)) & mask_4;

    // if (binaryToString(bc4, len_4) == "CGCC" || binaryToString(bc3, len_3)=="AAGA"){
    //   std::cout << binaryToString(bc4, len_4) << "\t" << binaryToString(bc3, len_3) << "\n";
    // }

    correct_34[bc4].second.add(bc3);
    correct_34[bc3].first.add(bc4);
  }

  std::streambuf *buf = nullptr;
  std::ofstream busf_out;

  if (!opt.stream_out)
  {
    busf_out.open(opt.output, std::ios::out | std::ios::binary);
    buf = busf_out.rdbuf();
  }
  else
  {
    buf = std::cout.rdbuf();
  }
  std::ostream bus_out(buf);

  bool outheader_written = false;

  nr = 0;
  BUSData bd;
  for (const auto &infn : opt.files)
  {
    std::streambuf *inbuf;
    std::ifstream inf;
    if (!opt.stream_in)
    {
      inf.open(infn.c_str(), std::ios::binary);
      inbuf = inf.rdbuf();
    }
    else
    {
      inbuf = std::cin.rdbuf();
    }
    std::istream in(inbuf);
    parseHeader(in, h);

    if (!outheader_written)
    {
      writeHeader(bus_out, h);
      outheader_written = true;
    }

    if (bclen == 0)
    {
      bclen = h.bclen;

      if (bclen != len_12 + len_34)
      {
        std::cerr << "Error: barcode length and whitelist length differ, barcodes = " << bclen << ", whitelist = " << len_12 + len_34 << std::endl
                  << "       check that your whitelist matches the technology used" << std::endl;

        exit(1);
      }
    }
    if (umilen == 0)
    {
      umilen = h.umilen;
    }

    int rc = 0;
    while (true)
    {
      in.read((char *)p, N * sizeof(BUSData));
      size_t rc = in.gcount() / sizeof(BUSData);
      if (rc == 0)
      {
        break;
      }
      nr += rc;

      for (size_t i = 0; i < rc; i++)
      {

        int n_12;
        int n_34;

        bd = p[i];

        uint64_t b = bd.barcode;
        uint64_t bc12 = b & mask_12;
        uint64_t bc34 = (b >> (2 * len_12)) & mask_34;

        auto it_12 = wbc_12.find(bc12);
        auto it_34 = wbc_34.find(bc34);

        it_12 != wbc_12.end() ? n_12 = 1 : n_12 = 0;
        it_34 != wbc_34.end() ? n_34 = 1 : n_34 = 0;

        if (n_12 == 1 && n_34 == 1)
        {
          stat_white++;
          bus_out.write((char *)&bd, sizeof(bd));
        }
        else if (n_12 + n_34 == 0 or n_12 + n_34 == 1) // either 12 or 34 == 0 or both ==0
        {                                              // TODO: there is code redundancy in this else statement, fix.
                                                       // CGCCAAGA GAGGGCAT top most
                                                       // CGCCAAGA AAGGGCAT next top most
                                                       // if (binaryToString(b, bclen) == "CGCCAAGAAAGGGCAT" || binaryToString(b, bclen) == "CGCCAAGAGAGGGCAT"){
                                                       //   std::cout << binaryToString(b, bclen) << "\t" << n_34 << "\t" << n_12 << "\n";
                                                       // }
          uint64_t bc12_corrected = bc12;
          uint64_t bc34_corrected = bc34;

          bool corrected_12_flag = false;
          bool corrected_34_flag = false;

          if (n_12 == 0)
          {
            uint64_t bc1 = b & mask_1;
            uint64_t bc2 = (b >> (2 * len_1)) & mask_2;

            uint64_t bc1_c = 0;
            uint64_t bc2_c = 0;

            int correct_1 = search_for_mismatch(correct_12[bc2].second, len_1, bc1, bc1_c);
            int correct_2 = search_for_mismatch(correct_12[bc1].first, len_2, bc2, bc2_c);

            int check_12 = correct_1 + correct_2;

            if (check_12 == 1)
            {
              n_12 = 1;
              corrected_12_flag = true;
              if (correct_1 == 1)
              {
                bc12_corrected = (bc2 << (2 * len_1)) | bc1_c; // 1 has been corrected, 2 is the same
              }
              else if (correct_2 == 1)
              {
                bc12_corrected = (bc2_c << (2 * len_1)) | bc1; // 2 has been corrected, 1 is the same
              }
            }
          }

          if (n_34 == 0)
          {
            uint64_t bc3 = bc34 & mask_3;
            uint64_t bc4 = (bc34 >> (2 * len_3)) & mask_4;

            uint64_t bc3_c = 0;
            uint64_t bc4_c = 0;

            int correct_3 = search_for_mismatch(correct_34[bc4].second, len_3, bc3, bc3_c);
            int correct_4 = search_for_mismatch(correct_34[bc3].first, len_4, bc4, bc4_c);

            int check_34 = correct_3 + correct_4;

            if (check_34 == 1)
            {
              n_34 = 1;
              corrected_34_flag = true;
              if (correct_3 == 1)
              {
                bc34_corrected = (bc4 << (2 * len_3)) | bc3_c; // 3 has been corrected, 4 is the same
              }
              else if (correct_4 == 1)
              {
                bc34_corrected = (bc4_c << (2 * len_3)) | bc3; // 4 has been corrected, 3 is the same
              }
            }
          }

          if (n_12 == 1 && n_34 == 1)
          {
            uint64_t b_corrected = bc34_corrected << (2 * len_12) | bc12_corrected;

            if (dump_bool)
            {
              if (bd.barcode != old_barcode)
              {
                of << binaryToString(bd.barcode, bclen) << "\t" << binaryToString(b_corrected, bclen) << "\n";
                old_barcode = bd.barcode;
              }
            }

            bd.barcode = b_corrected;
            bus_out.write((char *)&bd, sizeof(bd));

            if (corrected_12_flag && corrected_34_flag)
            { // if both were corrected once
              stat_corr_2++;
            }
            else if (corrected_12_flag || corrected_34_flag)
            { // if only one was corrected once
              stat_corr++;
            }
          }
          else
          {
            stat_uncorr++;
          }

          // AGCCAAGA GAGGGCAT
        }
      }
    }
  }

  std::cerr << "Processed " << nr << " BUS records" << std::endl
            << "In whitelist = " << stat_white << std::endl
            << "Corrected 1  = " << stat_corr << std::endl
            << "Corrected 2  = " << stat_corr_2 << std::endl
            << "Uncorrected  = " << stat_uncorr << std::endl;

  if (!opt.stream_out)
  {
    busf_out.close();
  }
  if (opt.dump_bool)
  {
    of.close(); // if of is open
  }

  delete[] p;
  p = nullptr;
}

void bustools_correct(Bustools_opt &opt)
{
  uint32_t bclen = 0;
  uint32_t wc_bclen = 0;
  uint32_t umilen = 0;
  BUSHeader h;
  size_t nr = 0;
  size_t N = 100000;
  BUSData *p = new BUSData[N];
  char magic[4];
  uint32_t version = 0;
  size_t stat_white = 0;
  size_t stat_corr = 0;
  size_t stat_uncorr = 0;
  uint64_t old_barcode;

  bool dump_bool = opt.dump_bool;

  std::ofstream of;
  if (dump_bool)
  {
    of.open(opt.dump);
  }

  std::ifstream wf(opt.whitelist, std::ios::in);
  std::string line;
  line.reserve(100);
  std::unordered_set<uint64_t> wbc;
  wbc.reserve(100000);
  uint32_t f = 0;
  while (std::getline(wf, line))
  {
    if (wc_bclen == 0)
    {
      wc_bclen = line.size();
    }
    uint64_t bc = stringToBinary(line, f);
    wbc.insert(bc);
  }
  wf.close();

  std::cerr << "Found " << wbc.size() << " barcodes in the whitelist" << std::endl;

  // split barcode into upper and lower half
  size_t bc2 = (wc_bclen + 1) / 2;

  std::vector<std::pair<Roaring, Roaring>> correct(1ULL << (2 * bc2)); // 4^(bc/2) possible barcodes

  uint64_t mask_size = (1ULL << (2 * bc2));
  uint64_t lower_mask = (1ULL << (2 * bc2)) - 1;
  uint64_t upper_mask = (1ULL << (2 * (wc_bclen - bc2))) - 1;
  for (uint64_t b : wbc)
  {
    uint64_t lb = b & lower_mask;
    uint64_t ub = (b >> (2 * bc2)) & upper_mask;

    correct[ub].second.add(lb);
    correct[lb].first.add(ub);
  }

  std::streambuf *buf = nullptr;
  std::ofstream busf_out;

  if (!opt.stream_out)
  {
    busf_out.open(opt.output, std::ios::out | std::ios::binary);
    buf = busf_out.rdbuf();
  }
  else
  {
    buf = std::cout.rdbuf();
  }
  std::ostream bus_out(buf);

  bool outheader_written = false;

  nr = 0;
  BUSData bd;
  for (const auto &infn : opt.files)
  {
    std::streambuf *inbuf;
    std::ifstream inf;
    if (!opt.stream_in)
    {
      inf.open(infn.c_str(), std::ios::binary);
      inbuf = inf.rdbuf();
    }
    else
    {
      inbuf = std::cin.rdbuf();
    }
    std::istream in(inbuf);
    parseHeader(in, h);

    if (!outheader_written)
    {
      writeHeader(bus_out, h);
      outheader_written = true;
    }

    if (bclen == 0)
    {
      bclen = h.bclen;

      if (bclen != wc_bclen)
      {
        std::cerr << "Error: barcode length and whitelist length differ, barcodes = " << bclen << ", whitelist = " << wc_bclen << std::endl
                  << "       check that your whitelist matches the technology used" << std::endl;

        exit(1);
      }
    }
    if (umilen == 0)
    {
      umilen = h.umilen;
    }

    int rc = 0;
    while (true)
    {
      in.read((char *)p, N * sizeof(BUSData));
      size_t rc = in.gcount() / sizeof(BUSData);
      if (rc == 0)
      {
        break;
      }
      nr += rc;

      for (size_t i = 0; i < rc; i++)
      {
        bd = p[i];
        auto it = wbc.find(bd.barcode);
        if (it != wbc.end())
        {
          stat_white++;
          bus_out.write((char *)&bd, sizeof(bd));
        }
        else
        {
          uint64_t b = bd.barcode;
          uint64_t lb = b & lower_mask;
          uint64_t ub = (b >> (2 * bc2)) & upper_mask;
          uint64_t lbc = 0, ubc = 0;
          int correct_lower = search_for_mismatch(correct[ub].second, bc2, lb, lbc);
          int correct_upper = search_for_mismatch(correct[lb].first, wc_bclen - bc2, ub, ubc);
          int nc = correct_lower + correct_upper;
          if (nc != 1)
          {
            stat_uncorr++;
          }
          else if (nc == 1)
          {
            if (correct_lower == 1)
            {
              uint64_t b_corrected = (ub << (2 * bc2)) | lbc;
              if (dump_bool)
              {
                if (bd.barcode != old_barcode)
                {
                  of << binaryToString(bd.barcode, bclen) << "\t" << binaryToString(b_corrected, bclen) << "\n";
                  old_barcode = bd.barcode;
                }
              }

              bd.barcode = b_corrected;
              bus_out.write((char *)&bd, sizeof(bd));
              stat_corr++;
            }
            else if (correct_upper == 1)
            {
              uint64_t b_corrected = (ubc << (2 * bc2)) | lb;
              if (dump_bool)
              {
                if (bd.barcode != old_barcode)
                {
                  of << binaryToString(bd.barcode, bclen) << "\t" << binaryToString(b_corrected, bclen) << "\n";
                  old_barcode = bd.barcode;
                }
              }

              bd.barcode = b_corrected;
              bus_out.write((char *)&bd, sizeof(bd));
              stat_corr++;
            }
          }
        }
      }
    }
  }

  std::cerr << "Processed " << nr << " BUS records" << std::endl
            << "In whitelist = " << stat_white << std::endl
            << "Corrected    = " << stat_corr << std::endl
            << "Uncorrected  = " << stat_uncorr << std::endl;

  if (!opt.stream_out)
  {
    busf_out.close();
  }

  if (opt.dump_bool)
  {
    of.close(); // if of is open
  }

  delete[] p;
  p = nullptr;
}