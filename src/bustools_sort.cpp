#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <queue>
#include <functional>
#include <thread>

#include "Common.hpp"
#include "BUSData.h"
#include "bustools_sort.h"

#include <time.h>

#define TP std::pair<BUSData, int>

//This code is for automatically creating the tmp directory supplied if it doesn't exist
#if defined(__MINGW32__) || defined(_MSC_VER)
//#include <filesystem> //once filesystem is acceptable for minGW, switch to that
#include "windows.h" //Needed for CreateDirectory


  void EnsureWindowsTempDirectoryExists(const Bustools_opt &opt) {
    //Make sure to create the tmp directory if it doesn't exist - writing temporary files fails otherwise in Windows
    //First get the directory - in theory, opt.temp_files can look like "tmp/x_" or just "x_" (or even nothing)
    //so we should find the last slash and make sure that directory exists

    std::size_t ind = opt.temp_files.rfind('/');
    std::size_t ind2 = opt.temp_files.rfind('\\');
    if (ind == std::string::npos)
    {
    ind = ind2;
    }
    else if (ind2 != std::string::npos)
    {
    //both valid, take the largest value (representing the last slash)
    ind = std::max(ind, ind2);
    }
    if (ind != std::string::npos)
    {
      auto dirName = opt.temp_files.substr(0, ind);
    //When our MinGW builds support c++17, change to std::filesystem
  
      //std::filesystem::path filepath = dirName;
      //if (!std::filesystem::is_directory(filepath)) 
      //{
      //    std::filesystem::create_directory(filepath);
      //}
    CreateDirectory(dirName.c_str(), NULL); //This will do nothing if the directory exists already
    }
  }

    //There is a bug in Windows, where bustools sort fails. The problem is that 
  //gcount for some reason fails here if too much is read and returns 0, even though
  //it succeeds. Could perhaps be a 32 bit issue somewhere, does size_t become 32 bits?
  //Anyway, this is a workaround that fixes the issue - does the same as the flag -m 100000000.
  //An interesting observation is that opt.max_memory is set to 1 << 32, which will become exactly
  //zero if truncated to 32 bits...
  size_t WindowsMaxMemory(size_t mem) {

    const size_t win_mem_max = 1e8;
    if (mem > win_mem_max) {
    mem = win_mem_max;
    }
    return mem;
  }
#else
  void EnsureWindowsTempDirectoryExists(const Bustools_opt &opt) {}
  size_t WindowsMaxMemory(size_t mem) {return mem;}

#endif





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
    if (a.pad == b.pad)
    {
      if (a.ec == b.ec)
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
        return a.ec < b.ec;
      }
    }
    else
    {
      return a.pad < b.pad;
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
    if (a.first.pad == b.first.pad)
    {
      if (a.first.ec == b.first.ec)
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
        return a.first.ec > b.first.ec;
      }
    }
    else
    {
      return a.first.pad > b.first.pad;
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

inline bool cmp5(const BUSData &a, const BUSData &b)
{
  if (a.flags == b.flags)
  {
    if (a.pad == b.pad)
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
      return a.pad < b.pad;
    }
  }
  else
  {
    return a.flags < b.flags;
  }
};

inline bool ncmp5(const TP &a, const TP &b)
{
  if (a.first.flags == b.first.flags)
  {
    if (a.first.pad == b.first.pad)
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
      return a.first.pad > b.first.pad;
    }
  }
  else
  {
    return a.first.flags > b.first.flags;
  }
};


void sort_bus_array(BUSData* busdata, size_t N, const int t, bool (*cmp)(const BUSData &, const BUSData &)) {
  //std::sort(busdata, busdata + N, cmp);
  if (t > 1 && N > 100000) {
    const size_t s = 256;
    std::vector<BUSData> samples, pivots;

    // samples = drawn from 0, s, 2s, ... , t*s
    samples.reserve(s*t);
    for (int i = 0; i < s*t; ++i) {
      samples.push_back(busdata[i * (N / (s*t))]);
    }
    std::sort(samples.begin(), samples.end(), cmp);

    pivots.reserve(t-1);
    // piviots are samples s, 2s, ... , (t-1)*s
    for (int i = 1; i < t; ++i) {
      pivots.push_back(samples[i * s]);
      
      //std::cerr << "pivot " << i << " = " << binaryToString(pivots[i-1].barcode, 16) << std::endl;
    }

    // buckets are locations of pivots after partitioning
    // partition i is between buckets[i] and buckets[i+1]
    std::vector<size_t> buckets(t+1, 0);
    buckets[0] = 0;
    buckets[t] = N;


    double partition_time = 0;
    clock_t start, end;
    start = clock();
    /*
    for (int i = 0; i < t-1; i++) {
      BUSData p = pivots[i];
      //std::cerr << "partitioning around " << binaryToString(p.barcode, 16) << std::endl;
      auto mid =  std::partition(busdata + buckets[i], busdata + N, [&p, &cmp](const BUSData &a) { return cmp(a, p); }) - busdata;
      buckets[i+1] = mid;
      //std::cerr << "bucket " << i << " has " << buckets[i+1] - buckets[i] << " elements, mid =  " << mid << std::endl;
    }
    //std::cerr << "bucket " << t-1 << " has " << buckets[t] - buckets[t-1] << " elements" << std::endl;
    */

    // partition the busdata based on the middle pivot
    std::function<void(int,int)> mid_partition = [&](int i, int j) { 
      if (j-i <= 1) {
        return;
      }
      size_t k = i + (j-i)/2;
      //std::cerr << "partitioning " << i << " to " << j << " with middle " << k << std::endl;
      //std::cerr << "buckets i and j are " << buckets[i] << " and " << buckets[j] << std::endl;
      BUSData p = pivots[k-1];
      //std::cerr << "pivot element is " << binaryToString(p.barcode, 16) << std::endl;
      buckets[k] = std::partition(busdata + buckets[i], busdata + buckets[j], [&p, &cmp](const BUSData &a) { return cmp(a, p); }) - busdata;
      //std::cerr << "bucket " << k << " is " <<  buckets[k] << std::endl;
      mid_partition(i, k);
      mid_partition(k, j);
    };

    mid_partition(0, t);
    

    /*
    // verify that the pivots are sorted
    for (int i = 0; i < t-2; i++) {
      if (!cmp(pivots[i], pivots[i+1])) {
        std::cerr << "pivot " << i << " is not smaller than pivot " << i+1 << std::endl;
        exit(1);
      }
    }
    
    for (int i = 0; i < t; i++) {
      std::cerr << "bucket " << i << " at " << buckets[i] <<  " has " << buckets[i+1] - buckets[i] << " elements" << std::endl;
    }
    
    //verify that each partition is smaller than the pivot
    for (int i = 0; i < t; i++) {
      for (size_t j = buckets[i]; j < buckets[i+1]; j++) {
        if (i < t-1 && !cmp(busdata[j], pivots[i])) {
          std::cerr << "partition " << i << " has an element larger than the pivot" << std::endl;
          std::cerr << "element " << j << " = " << binaryToString(busdata[j].barcode, 16) << std::endl;
          std::cerr << "pivot " << i << " = " << binaryToString(pivots[i].barcode, 16) << std::endl;

          exit(1);
        }
        if (i > 0 && cmp(busdata[j], pivots[i-1])) {
          std::cerr << "partition " << i << " has an element smaller than the next pivot" << std::endl;
          std::cerr << "element " << j << " = " << binaryToString(busdata[j].barcode, 16) << std::endl;
          std::cerr << "pivot " << i-1 << " = " << binaryToString(pivots[i-1].barcode, 16) << std::endl;
          exit(1);
        }
      }
    }
    */

    
    end = clock();
    partition_time  += ((double) (end - start)) / CLOCKS_PER_SEC;
    std::cerr << "partition time: " << partition_time << "s" << std::endl;



    // sort each bucket
    std::vector<std::thread> workers;
    for (int i = 0; i < t; ++i) {
      workers.push_back(std::thread([&busdata, &buckets, &cmp, i]() {
        //std::cerr << "sorting bucket " << i << " with " << buckets[i] <<  " to " << buckets[i+1]<< std::endl;
        std::sort(busdata + buckets[i], busdata + buckets[i+1], cmp);
      }));
    }

    for (auto &w : workers) {
      w.join();
    }


  } else {
    std::sort(busdata, busdata + N, cmp);
  }

}

void bustools_sort(const Bustools_opt &opt)
{
  auto mem = WindowsMaxMemory(opt.max_memory);

  BUSHeader h;
  size_t N = mem / sizeof(BUSData);
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
  case SORT_F_BC:
    cmp = &cmp5;
    ncmp = &ncmp5;
    break;
  default:
    std::cerr << "ERROR: Unknown sort type" << std::endl;
    exit(1);
  }


  

  size_t sc = 0; // number of records read
  double sorting_time = 0;
  int tmp_file_no = 0;

  // only use a single buffer if we are reading from stdin or if we have a single file
  bool all_in_buffer = opt.stream_in || opt.files.size() == 1;


  const auto collapse_and_write = [&](BUSData *p, size_t rc, std::ostream &outf) {
    for (size_t i = 0; i < rc;) {
      size_t j = i + 1;
      uint32_t c = p[i].count;
      auto ec = p[i].ec;
      for (; j < rc; j++) {
        if (p[i].barcode == p[j].barcode && p[i].UMI == p[j].UMI && p[i].ec == p[j].ec && p[i].flags == p[j].flags && p[i].pad == p[j].pad) {
          c += p[j].count;
        } else {
          break;
        }
      }
      // merge identical things
      p[i].count = c;
      outf.write((char *)(&(p[i])), sizeof(p[i]));
      // increment
      i = j;
    }
  };

  // open the correct output stream
  std::ofstream of;
  std::streambuf *buf = nullptr;
  if (!opt.stream_out) {
      of.open(opt.output, std::ios::out | std::ios::binary);
      buf = of.rdbuf();
  } else {
      buf = std::cout.rdbuf();
  }
  std::ostream busf_out(buf);

  // measure time spent reading input
  clock_t start,end;
  double reading_time = 0;
  double writing_time = 0;

  for (const auto &infn : opt.files) {
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

    int rc = 1;


    while (in.good()) {
      
      start = clock();
      in.read((char *)p, N * sizeof(BUSData));
      size_t rc = in.gcount() / sizeof(BUSData);
      end = clock();
      reading_time += ((double) (end - start)) / CLOCKS_PER_SEC;

      // no records read, we are done
      if (rc == 0) {
        break;
      }

      // records did not fit in buffer
      if (rc >= N) {
        all_in_buffer = false;
      }
      
      // now sort the data
      start = clock();
      //std::sort(p, p + rc, cmp);
      sort_bus_array(p, rc, opt.threads, cmp);
      end = clock();
      sorting_time += ((double) (end - start)) / CLOCKS_PER_SEC;

      sc += rc;

      if (all_in_buffer) {
        std::cerr << " all fits in buffer" << std::endl;
        // single file or stream, all data fits in buffer, write directly to output
        start = clock();
        writeHeader(busf_out, h);
        collapse_and_write(p, rc, busf_out);
        end = clock();
        writing_time = ((double) (end - start)) / CLOCKS_PER_SEC;
      } else {
        // need to sort in chunks
        // write the output
        std::ofstream outf(opt.temp_files + std::to_string(tmp_file_no), std::ios::binary);
        writeHeader(outf, h);

        collapse_and_write(p, rc, outf);
    
        outf.close();
        tmp_file_no++;
      }
      
    }
  }
  delete[] p;
  p = nullptr;

  std::cerr << "Read in " << sc << " BUS records" << std::endl;
  
  std::cerr << "reading time " << reading_time << "s" << std::endl;
  std::cerr << "sorting time " << sorting_time << "s" << std::endl;
  std::cerr << "writing time " << writing_time << "s" << std::endl;



  
  if (!all_in_buffer) {
    writeHeader(busf_out, h);

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
