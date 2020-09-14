#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <queue>
#include <functional>

#include "Common.hpp"
#include "BUSData.h"
#include "bustools_sort.h"

#define TP std::pair<BUSData, int>

inline bool cmp1(const BUSData &a, const BUSData &b)
{
  if (a.barcode == b.barcode)
  {
    if (a.UMI == b.UMI)
    {
      if (a.ec == b.ec)
      {
        return a.flags < b.flags;
      }
      else
      {
        return a.ec < b.ec;
      }
    }
    else
    {
      return a.UMI < b.UMI;
    }
  }
  else
  {
    return a.barcode < b.barcode;
  }
};

inline bool ncmp1(const TP &a, const TP &b)
{
  if (a.first.barcode == b.first.barcode)
  {
    if (a.first.UMI == b.first.UMI)
    {
      if (a.first.ec == b.first.ec)
      {
        if (a.first.flags == b.first.flags)
        {
          return a.second > b.second;
        }
        else
        {
          return a.first.flags > b.first.flags;
        }
      }
      else
      {
        return a.first.ec > b.first.ec;
      }
    }
    else
    {
      return a.first.UMI > b.first.UMI;
    }
  }
  else
  {
    return a.first.barcode > b.first.barcode;
  }
};

inline bool cmp2(const BUSData &a, const BUSData &b)
{
  if (a.UMI == b.UMI)
  {
    if (a.barcode == b.barcode)
    {
      return a.ec < b.ec;
    }
    else
    {
      return a.barcode < b.barcode;
    }
  }
  else
  {
    return a.UMI < b.UMI;
  }
};

inline bool ncmp2(const TP &a, const TP &b)
{
  if (a.first.UMI == b.first.UMI)
  {
    if (a.first.barcode == b.first.barcode)
    {
      if (a.first.ec == b.first.ec)
      {
        return a.second > b.second;
      }
      else
      {
        return a.first.ec > b.first.ec;
      }
    }
    else
    {
      return a.first.barcode > b.first.barcode;
    }
  }
  else
  {
    return a.first.UMI > b.first.UMI;
  }
};

inline bool cmp3(const BUSData &a, const BUSData &b)
{
  if (a.flags == b.flags)
  {
    if (a.ec == b.ec)
    {
      if (a.pad == b.pad)
      {
        if (a.barcode == b.barcode)
        {
          return a.UMI < b.UMI;
        }
        else
        {
          return a.barcode < b.barcode;
        }
      }
      else
      {
        return a.pad < b.pad;
      }
    }
    else
    {
      return a.ec < b.ec;
    }
  }
  else
  {
    return a.flags < b.flags;
  }
};

inline bool ncmp3(const TP &a, const TP &b)
{
  if (a.first.flags == b.first.flags)
  {
    if (a.first.ec == b.first.ec)
    {
      if (a.first.pad == b.first.pad)
      {
        if (a.first.barcode == b.first.barcode)
        {
          if (a.first.UMI == b.first.UMI)
          {
            return a.second > b.second;
          }
          else
          {
            return a.first.UMI > b.first.UMI;
          }
        }
        else
        {
          return a.first.barcode > b.first.barcode;
        }
      }
      else
      {
        return a.first.pad > b.first.pad;
      }
    }
    else
    {
      return a.first.ec > b.first.ec;
    }
  }
  else
  {
    return a.first.flags > b.first.flags;
  }
};

inline bool cmp4(const BUSData &a, const BUSData &b)
{
  if (a.count == b.count)
  {
    if (a.barcode == b.barcode)
    {
      if (a.UMI == b.UMI)
      {
        return a.ec < b.ec;
      }
      else
      {
        return a.UMI < b.UMI;
      }
    }
    else
    {
      return a.barcode < b.barcode;
    }
  }
  else
  {
    return a.count < b.count;
  }
};

inline bool ncmp4(const TP &a, const TP &b)
{
  if (a.first.count == b.first.count)
  {
    if (a.first.barcode == b.first.barcode)
    {
      if (a.first.UMI == b.first.UMI)
      {
        if (a.first.ec == b.first.ec)
        {
          return a.second > b.second;
        }
        else
        {
          return a.first.ec > b.first.ec;
        }
      }
      else
      {
        return a.first.UMI > b.first.UMI;
      }
    }
    else
    {
      return a.first.barcode > b.first.barcode;
    }
  }
  else
  {
    return a.first.count > b.first.count;
  }
};

void bustools_sort(const Bustools_opt &opt)
{
  BUSHeader h;
  size_t N = opt.max_memory / sizeof(BUSData);
  BUSData *p = new BUSData[N];
  char magic[4];
  uint32_t version = 0;

  int no_temp_files = 0;

  bool (*cmp)(const BUSData &, const BUSData &);
  bool (*ncmp)(const TP &a, const TP &b);
  switch (opt.type)
  {
  case SORT_BC:
    cmp = &cmp1;
    ncmp = &ncmp1;
    break;
  case SORT_UMI:
    cmp = &cmp2;
    ncmp = &ncmp2;
    break;
  case SORT_F:
    cmp = &cmp3;
    ncmp = &ncmp3;
    break;
  case SORT_COUNT:
    cmp = &cmp4;
    ncmp = &ncmp4;
    break;
  default:
    std::cerr << "ERROR: Unknown sort type" << std::endl;
    exit(1);
  }

  size_t sc = 0;
  int tmp_file_no = 0;
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

    int rc = 1;

    while (in.good())
    {
      // read as much as we can
      in.read((char *)p, N * sizeof(BUSData));
      size_t rc = in.gcount() / sizeof(BUSData);
      if (rc == 0)
      {
        break;
      }
      // now sort the data
      std::sort(p, p + rc, cmp);
      sc += rc;

      // write the output
      std::ofstream outf(opt.temp_files + std::to_string(tmp_file_no), std::ios::binary);
      writeHeader(outf, h);

      for (size_t i = 0; i < rc;)
      {
        size_t j = i + 1;
        uint32_t c = p[i].count;
        auto ec = p[i].ec;
        for (; j < rc; j++)
        {
          if (p[i].barcode == p[j].barcode && p[i].UMI == p[j].UMI && p[i].ec == p[j].ec && p[i].flags == p[j].flags && p[i].pad == p[j].pad)
          {
            c += p[j].count;
          }
          else
          {
            break;
          }
        }
        // merge identical things
        p[i].count = c;
        outf.write((char *)(&(p[i])), sizeof(p[i]));
        // increment
        i = j;
      }

      outf.close();
      tmp_file_no++;
    }
  }
  delete[] p;
  p = nullptr;

  std::cerr << "Read in " << sc << " BUS records" << std::endl;

  std::streambuf *buf = nullptr;
  std::ofstream of;

  if (!opt.stream_out)
  {
    of.open(opt.output, std::ios::out | std::ios::binary);
    buf = of.rdbuf();
  }
  else
  {
    buf = std::cout.rdbuf();
  }
  std::ostream busf_out(buf);

  writeHeader(busf_out, h);

  // todo: skip writing to disk if it fits in memory
  if (tmp_file_no == 1)
  {
    size_t M = N / 8;
    p = new BUSData[M];
    std::ifstream in(opt.temp_files + "0", std::ios::binary);
    BUSHeader tmp;
    parseHeader(in, tmp);
    while (in.good())
    {
      // read as much as we can
      in.read((char *)p, M * sizeof(BUSData));
      size_t rc = in.gcount() / sizeof(BUSData);
      if (rc == 0)
      {
        break;
      }
      busf_out.write((char *)p, rc * sizeof(BUSData));
    }
    in.close();
    std::remove((opt.temp_files + "0").c_str());
  }
  else
  {
    // TODO: test if replacing with k-way merge is better
    // adapted from https://github.com/arq5x/kway-mergesort/blob/master/kwaymergesort.h
    int k = tmp_file_no;
    size_t M = N / (k);
    //std::memset(p, 0, N*sizeof(BUSData));
    std::vector<std::ifstream> bf(k);
    for (int i = 0; i < k; i++)
    {
      bf[i].open((opt.temp_files + std::to_string(i)).c_str(), std::ios::binary);
      BUSHeader tmp;
      parseHeader(bf[i], tmp);
    }

    std::priority_queue<TP, std::vector<TP>, std::function<bool(const TP &a, const TP &b)>> pq(ncmp);
    BUSData t;
    for (int i = 0; i < k; i++)
    {
      bf[i].read((char *)&t, sizeof(t));
      pq.push({t, i});
    }

    BUSData curr = pq.top().first;
    curr.count = 0; // we'll count this again in the first loop
    while (!pq.empty())
    {
      TP min = pq.top();

      pq.pop();
      // process the data
      BUSData &m = min.first;
      int i = min.second;
      if (m.barcode == curr.barcode && m.UMI == curr.UMI && m.ec == curr.ec && m.flags == curr.flags && m.pad == curr.pad)
      {
        // same data, increase count
        curr.count += m.count;
      }
      else
      {

        // new data let's output curr, new curr is m
        if (curr.count != 0)
        {
          busf_out.write((char *)&curr, sizeof(curr));
        }
        curr = m;
      }
      // read next from stream
      if (bf[i].good())
      {
        bf[i].read((char *)&t, sizeof(t));
        if (bf[i].gcount() > 0)
        {
          pq.push({t, i});
        }
      }
    }

    if (curr.count > 0)
    {
      // write out remaining straggler
      busf_out.write((char *)&curr, sizeof(curr));
    }

    // remove intermediary files
    for (int i = 0; i < k; i++)
    {
      bf[i].close();
      std::remove((opt.temp_files + std::to_string(i)).c_str());
    }
  }

  if (!opt.stream_out)
  {
    of.close();
  }
}

void bustools_sort_orig(const Bustools_opt &opt)
{
  BUSHeader h;
  std::vector<BUSData> b;
  size_t N = 100000;
  BUSData *p = new BUSData[N];
  char magic[4];
  uint32_t version = 0;
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

    int rc = 1;
    while (true)
    {
      in.read((char *)p, N * sizeof(BUSData));
      size_t rc = in.gcount() / sizeof(BUSData);
      if (rc == 0)
      {
        break;
      }
      // todo, reserve max memory
      b.insert(b.end(), p, p + rc);
    }
  }

  delete[] p;
  p = nullptr;
  std::cerr << "Read in " << b.size() << " BUS records" << std::endl;

  // todo: replace with radix sort
  std::sort(b.begin(), b.end(), [&](const BUSData &a, const BUSData &b) {
                                    if (a.barcode == b.barcode) {
                                      if (a.UMI == b.UMI) {
                                        return a.ec < b.ec;
                                      } else {
                                        return a.UMI < b.UMI;
                                      } 
                                    } else { 
                                      return a.barcode < b.barcode;
                                    } });
  std::cerr << "All sorted" << std::endl;

  std::streambuf *buf = nullptr;
  std::ofstream of;

  if (!opt.stream_out)
  {
    of.open(opt.output, std::ios::out | std::ios::binary);
    buf = of.rdbuf();
  }
  else
  {
    buf = std::cout.rdbuf();
  }
  std::ostream busf_out(buf);

  writeHeader(busf_out, h);

  size_t n = b.size();
  for (size_t i = 0; i < n;)
  {
    size_t j = i + 1;
    uint32_t c = b[i].count;
    auto ec = b[i].ec;
    for (; j < n; j++)
    {
      if (b[i].barcode != b[j].barcode || b[i].UMI != b[j].UMI || b[i].ec != b[j].ec)
      {
        break;
      }
      c += b[j].count;
    }
    // merge identical things
    b[i].count = c;
    busf_out.write((char *)(&(b[i])), sizeof(b[i]));
    // increment
    i = j;
  }

  if (!opt.stream_out)
  {
    of.close();
  }
}
