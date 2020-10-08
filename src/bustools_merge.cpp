#include <iostream>
#include <fstream>
#include <algorithm>
#include <queue>
#include <functional>
#include <unordered_map>
#include <unordered_set>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_merge.h"

#define BPair std::pair<BUSData, BUSData>
#define BP std::pair<BPair, int> // bus pair, file that they came from
#define BUSRange std::pair<BUSData, std::pair<int32_t, int32_t>>
#define TP std::pair<BUSData, int>
#define Range std::pair<uint32_t, uint32_t>

inline std::vector<int32_t> intersect_vecs(const std::vector<int32_t> &x, const std::vector<int32_t> &y)
{
  std::vector<int32_t> v;
  auto a = x.begin();
  auto b = y.begin();
  while (a != x.end() && b != y.end())
  {
    if (*a < *b)
    {
      ++a;
    }
    else if (*b < *a)
    {
      ++b;
    }
    else
    {
      v.push_back(*a);
      ++a;
      ++b;
    }
  }
  return v;
}

inline std::vector<int32_t> get_tids(const BUSHeader &oh, const std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> &ecmapinv, const int32_t &eid)
{

  std::vector<int32_t> tids = oh.ecs[eid];

  // std::sort(tids.begin(), tids.end());
  // tids.erase(std::unique(tids.begin(), tids.end()), tids.end());

  return tids;
}

inline void print_bd(const BUSData &bd, const size_t bclen, const size_t umilen)
{
  std::cout << binaryToString(bd.barcode, bclen) << "\t" << binaryToString(bd.UMI, umilen) << "\t" << bd.ec << "\t" << bd.count << "\t" << bd.flags << "\t" << bd.pad << std::endl;
}

// vocab
// txn: string representing transcript str
// tid: number indexing the transcript in the file
// ecs: the set of tids (set of transcript indexes)
// eid: the number indexing the equivalence class
// keep in mind this could be wrt original file or
// merged file

void bustools_merge_different_index(const Bustools_opt &opt)
{
  // read in transcripts.txt
  std::ifstream ifn(opt.count_txp);
  std::string txn;
  int32_t tid;
  std::unordered_map<std::string, int32_t> txn_tid;
  std::vector<int32_t> tids;

  // insert tids into a vector
  while (ifn >> txn)
  {
    auto ok = txn_tid.insert({txn, tid});
    if (ok.second)
    {
      tids.push_back(tid);
      tid++;
    }
    else
    {
      tids.push_back(ok.first->second);
    }
  }
  ifn.close();
  std::cout << "[info] parsed transcripts.txt" << std::endl;

  // read in matrix.ec
  BUSHeader h, bh;
  parseECs(opt.count_ecs, h);
  // put the ecs into a ecmap inv
  std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> ecmapinv;

  for (int32_t ec; ec < h.ecs.size(); ec++)
  {
    ecmapinv.insert({h.ecs[ec], ec});
  }
  std::cout << "[info] parsed matrix.ec" << std::endl;

  // read in BUS Records grouped by flags
  std::streambuf *inbuf;
  std::ifstream inf;
  if (!opt.stream_in)
  {
    inf.open(opt.files[0].c_str(), std::ios::binary);
    inbuf = inf.rdbuf();
  }
  else
  {
    inbuf = std::cin.rdbuf();
  }

  std::istream in(inbuf);
  parseHeader(in, bh);
  uint32_t bclen = bh.bclen;
  uint32_t umilen = bh.umilen;
  int rc = 1;

  // implement std queue

  std::queue<BUSData> queue;
  BUSData bd, prev, curr;
  int N = 1024;
  size_t nr = 0, nw = 0;
  for (int i = 0; i < N; i++)
  {
    in.read((char *)&bd, sizeof(BUSData));
    queue.push(bd);
    nr++;
  }

  prev = queue.front();
  queue.pop();
  curr = queue.front();
  queue.pop();

  std::unordered_set<int32_t> c;
  std::vector<std::unordered_set<int32_t>> elem_sets;

  std::ofstream outf(opt.output);
  writeHeader(outf, bh);

  while (true)
  {
    // print_bd(prev, bclen, umilen);
    elem_sets.clear();
    c.clear();
    while (prev.flags == curr.flags && prev.barcode == curr.barcode && prev.UMI == curr.UMI)
    {
      // print_bd(curr, bclen, umilen);
      if (curr.ec != prev.ec)
      {
        if (c.size()) // if c has elements in the set
        {
          elem_sets.push_back(c);
        }
      }

      if (c.find(curr.ec) != c.end())
      {
        c.erase(curr.ec);
      }
      else
      {
        c.insert(curr.ec);
      }

      if (queue.empty())
      {
        break;
      }
      // read in more data
      prev = std::move(curr);
      curr = queue.front();
      queue.pop();
      if (in.good())
      {
        in.read((char *)&bd, sizeof(BUSData));
        queue.push(bd);
        nr++;
      }

      //print_bd(curr, bclen, umilen);
    }

    // do the intersection
    if (elem_sets.size())
    {

      //there is a weird segfault when printing this out..
      // std::cout << "before intersecting TIDS" << std::endl;
      // for (const std::unordered_set<int32_t> &e : elem_sets)
      // {
      //   if (!e.empty())
      //   {
      //     for (const int32_t &eid : e)
      //     {
      //       try
      //       {
      //         std::cout << h.ecs.size() << ": " << eid << ", ";
      //         if (eid > h.ecs.size())
      //         {
      //           std::cout << "--------SCREAM!!!!!!" << std::endl;
      //         }
      //         // std::vector<int32_t> tds = h.ecs[eid];
      //         // /* code */
      //         // for (auto &t : tds)
      //         // {
      //         //   std::cout << t << ", ";
      //         // }
      //       }
      //       catch (const std::exception &e)
      //       {
      //         std::cerr << e.what() << '\n';
      //         break;
      //       }
      //     }
      //   }
      //   std::cout << std::endl;
      // }

      std::vector<int32_t> tids_per_elem, prev_tids;
      bool single = true;
      std::vector<int32_t> tids;
      for (int32_t i = 0; i < elem_sets.size(); i++)
      {
        if (elem_sets[i].size())
        {
          for (const int32_t &e : elem_sets[i])
          {
            tids = h.ecs[e];
            tids_per_elem.insert(tids_per_elem.end(), tids.begin(), tids.end());
            tids.clear();
          }
          std::sort(tids_per_elem.begin(), tids_per_elem.end());
          tids_per_elem.erase(std::unique(tids_per_elem.begin(), tids_per_elem.end()), tids_per_elem.end());

          if (i == 0)
          {
            prev_tids = std::move(tids_per_elem);
          }
          else
          {
            prev_tids = intersect_vecs(prev_tids, tids_per_elem); // problem may be here..
            tids_per_elem.clear();
          }
          if (!prev_tids.size()) // an intermediary intersection is empty
          {
            break;
          }
        }
      }

      // std::cout << "After intersecting TIDS" << std::endl;
      // for (const int32_t &t : prev_tids)
      // {
      //   std::cout << t << ", ";
      // }
      // std::cout << std::endl;

      if (prev_tids.size()) // if there are tids left over write the record
      {
        // we find the ec associated with prev_tids
        auto findec = ecmapinv.find(prev_tids);
        // if the ec doesnt exist, make a new one
        if (findec == ecmapinv.end())
        {
          //std::cout << "end of ecmapinv" << std::endl;
          prev.ec = ecmapinv.size();
          h.ecs.push_back(prev_tids);
          ecmapinv.insert({prev_tids, prev.ec});
          // std::cout << "Making new ecs: " << prev.ec << ", " << ecmapinv.size() << ", " << h.ecs.size() << std::endl;
        }
        else // if it does exist, assign it
        {
          //std::cout << "found ecmap" << std::endl;
          prev.ec = findec->second;
        }
        // write the busrecord
        ///std::cout << "ec: " << BR.ec << std::endl;
        prev.count = 1;
        prev.pad = 0;
        outf.write((char *)&prev, sizeof(prev));
        nw++;
      }
      tids_per_elem.clear();
      prev_tids.clear();
    }

    // read in more data
    if (queue.empty())
    {
      break;
    }
    prev = std::move(curr);
    curr = queue.front();
    queue.pop();
    if (in.good())
    {
      in.read((char *)&bd, sizeof(BUSData));
      queue.push(bd);
      nr++;
    }
  }
  // std::cout << "end" << std::endl;
  std::cerr << "bus records read:    " << nr << std::endl;
  std::cerr << "bus records written: " << nw << std::endl;
}
