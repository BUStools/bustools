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
  size_t bc4 = (bc2 / 2);

  std::vector<size_t> split_bclen(4);

  std::vector<std::pair<Roaring, Roaring>> correct(1ULL << (2 * bc2)); // 4^(bc/2) possible barcodes

  std::vector<std::pair<Roaring, Roaring>> lower_correct(1ULL << (2 * bc4)); // 4^(bc/4) possible barcodes
  std::vector<std::pair<Roaring, Roaring>> upper_correct(1ULL << (2 * bc4)); // 4^(bc/4) possible barcodes

  uint64_t mask_size = (1ULL << (2 * bc2));
  uint64_t lower_mask = (1ULL << (2 * bc2)) - 1;
  uint64_t upper_mask = (1ULL << (2 * (wc_bclen - bc2))) - 1;

  uint64_t l_lower_mask = (1ULL << (2 * bc4)) - 1;
  uint64_t u_lower_mask = (1ULL << (2 * bc4)) - 1;

  uint64_t l_upper_mask = (1ULL << (2 * (bc2 - bc4))) - 1;
  uint64_t u_upper_mask = (1ULL << (2 * (bc2 - bc4))) - 1;

  // barcode is split into 4 peices of equal size, check each half in the following way:
  // AGCC AAGA GAGG GCAT
  // UU   LU    UL  LL
  // U: upper, L: lower
  // create all hamming 1 distance variants of the LL half and check to see if
  // any one of those variants (concatenated with the UL) is in the "half"-whitelist
  // Do the same for the UU and LU half.

  for (uint64_t b : wbc)
  {
    uint64_t lb = b & lower_mask;
    uint64_t ub = (b >> (2 * bc2)) & upper_mask;

    uint64_t llb = b & l_lower_mask;
    uint64_t ulb = (b >> (2 * bc4)) & u_lower_mask;

    uint64_t lub = (b >> (2 * (2 * bc4))) & l_upper_mask;
    uint64_t uub = (b >> (3 * (2 * bc4))) & u_upper_mask;

    // std::cout << "ub: " << binaryToString(ub, wc_bclen/2) << "\t" << "lb: " << binaryToString(lb, wc_bclen/2) << "\n";

    correct[ub].second.add(lb);
    correct[lb].first.add(ub);

    lower_correct[ulb].second.add(llb);
    lower_correct[llb].first.add(ulb);

    upper_correct[uub].second.add(lub);
    upper_correct[lub].first.add(uub);
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
        { // TODO: there is code redundancy in this else statement, fix.
          uint64_t b = bd.barcode;

          uint64_t lb = b & lower_mask;
          uint64_t ub = (b >> (2 * bc2)) & upper_mask;

          uint64_t upper_corrected = ub;
          uint64_t lower_corrected = lb;

          uint64_t llb = b & l_lower_mask;
          uint64_t ulb = (b >> (2 * bc4)) & u_lower_mask;

          uint64_t lub = (b >> (2 * (2 * bc4))) & l_upper_mask;
          uint64_t uub = (b >> (3 * (2 * bc4))) & u_upper_mask;

          uint64_t lbc = 0, ubc = 0;
          uint64_t llbc = 0, ulbc = 0;
          uint64_t lubc = 0, uubc = 0;

          // std::cout << binaryToString(uub, bc4) << "\t" << binaryToString(lub, bc4) << "\t"
          //           << binaryToString(ulb, bc4) << "\t" << binaryToString(llb, bc4) << "\n";

          int correct_l_lower = search_for_mismatch(lower_correct[ulb].second, bc4, llb, llbc);
          int correct_u_lower = search_for_mismatch(lower_correct[llb].first, bc4, ulb, ulbc);

          int nc_lower = correct_l_lower + correct_u_lower;

          if (nc_lower == 1)
          {
            if (correct_l_lower == 1)
            {
              lower_corrected = (ulb << (2 * bc4)) | llbc;
            }
            else if (correct_u_lower == 1)
            {
              lower_corrected = (ulbc << (2 * bc4)) | llb;
            }
          }

          int correct_l_upper = search_for_mismatch(upper_correct[uub].second, bc4, lub, lubc);
          int correct_u_upper = search_for_mismatch(upper_correct[lub].first, bc4, uub, uubc);

          int nc_upper = correct_l_upper + correct_u_upper;

          if (nc_upper == 1)
          {
            if (correct_l_upper == 1)
            {
              upper_corrected = (uub << (2 * bc4)) | lubc;
            }
            else if (correct_u_upper == 1)
            {
              upper_corrected = (uubc << (2 * bc4)) | lub;
            }
          }

          uint64_t b_corrected = upper_corrected << (2 * bc2) | lower_corrected;
          auto it = wbc.find(b_corrected);

          // int correct_lower = search_for_mismatch(correct[upper_corrected].second,bc2,lb,lbc);
          // int correct_upper = search_for_mismatch(correct[lower_corrected].first,wc_bclen - bc2,ub,ubc);
          // int nc = correct_lower + correct_upper;

          int nc = nc_lower + nc_upper;
          // AGCCAAGA GAGGGCAT

          if (nc == 0 || nc > 2)
          {
            stat_uncorr++;
          }
          else if (nc == 1)
          {
            if (dump_bool)
            {
              if (it != wbc.end())
              {
                bd.barcode = b_corrected;
                stat_corr++;
                bus_out.write((char *)&bd, sizeof(bd));

                if (bd.barcode != old_barcode)
                {
                  of << binaryToString(bd.barcode, bclen) << "\t" << binaryToString(b_corrected, bclen) << "\n";
                  old_barcode = bd.barcode;
                }
              }
            }
          }
            else if (nc == 2)
            {
              if (nc_lower==1 & nc_upper==1)
              {
                if (dump_bool)
                {
                  if (it != wbc.end())
                  {
                    bd.barcode = b_corrected;
                    stat_corr_2++;
                    bus_out.write((char *)&bd, sizeof(bd));

                    if (bd.barcode != old_barcode)
                    {
                      of << binaryToString(bd.barcode, bclen) << "\t" << binaryToString(b_corrected, bclen) << "\n";
                      old_barcode = bd.barcode;
                    }
                  }
                }
              }
            }
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
