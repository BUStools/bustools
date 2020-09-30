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

  std::sort(tids.begin(), tids.end());
  tids.erase(std::unique(tids.begin(), tids.end()), tids.end());

  return tids;
}

inline void print_bd(const BUSData &bd, const size_t bclen, const size_t umilen)
{
  std::cout << binaryToString(bd.barcode, bclen) << "\t" << binaryToString(bd.UMI, umilen) << "\t" << bd.ec << "\t" << bd.count << "\t" << bd.flags << "\t" << bd.pad << std::endl;
}

inline bool itv_cmp_left(const BUSRange &x, const BUSRange &y)
{
  if (x.second.first == y.second.first)
  {
    return x.second.second < y.second.second;
  }
  else
  {
    return x.second.first < y.second.first;
  }
}

inline bool itv_cmp_right(const BUSRange &x, BUSRange &y)
{
  if (x.second.second == y.second.second)
  {
    return x.second.first < y.second.first;
  }
  else
  {
    return x.second.second < y.second.second;
  }
}

// inline bool ncmp(const BP &a, const BP &b)
// {
//   if (a.first.flags == b.first.flags)
//   {
//     if (a.first.ec == b.first.ec)
//     {
//       return a.second > b.second;
//     }
//     else
//     {
//       return a.first.ec > b.first.ec;
//     }
//   }
//   else
//   {
//     return a.first.flags > b.first.flags;
//   }
// }

inline bool pqcmp(const BP &a, const BP &b)
{
  if (a.first.first.flags == b.first.first.flags)
  {
    return a.second > b.second;
  }
  else
  {
    return a.first.first.flags > b.first.first.flags;
  }
}

inline std::pair<std::vector<std::unordered_set<int32_t>>, std::vector<std::pair<int32_t, int32_t>>> find_elem(std::vector<BUSRange> &l, std::vector<BUSRange> &r)
{
  std::vector<std::unordered_set<int32_t>> d;
  std::vector<std::pair<int32_t, int32_t>> bounds;

  // initializers
  int lidx, ridx;
  int32_t start, stop;
  std::unordered_set<int32_t> c;
  // declarations
  lidx = 0;
  start = l[lidx].second.first;

  for (int ridx = 0; ridx < r.size(); ridx++)
  {

    c.insert(r[ridx].first.ec);

    while (lidx < l.size() && l[lidx].second.first <= r[ridx].second.second)
    {
      stop = l[lidx].second.first;

      c.insert(l[lidx].first.ec);

      if (stop > start)
      {
        c.erase(l[lidx].first.ec);
        d.push_back(c);
        bounds.push_back({start, stop});
        start = stop;
        c.insert(l[lidx].first.ec);
      }
      lidx += 1;
    }

    if (start < r[ridx].second.second)
    {
      stop = r[ridx].second.second;
      d.push_back(c);
      bounds.push_back({start, stop});
      start = stop;
    }
    c.erase(r[ridx].first.ec);
  }

  return {d, bounds};
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
  int nf = opt.files.size();

  std::vector<std::ifstream> bf(nf); // vector of bus file streams
  std::vector<BUSHeader> vh;         // vector of bus headers

  // populate the headers (we ignore ec for now)
  for (int i = 0; i < nf; i++)
  {
    bf[i].open((opt.files[i] + "/output.bus").c_str(), std::ios::binary);
    BUSHeader h;
    parseHeader(bf[i], h); // parse the header into h

    parseECs(opt.files[i] + "/matrix.ec", h); // parse ecs
    vh.push_back(std::move(h));               // place the ecs into h
  }
  std::cout << "[info] parsed output.bus files" << std::endl;

  // parse the transcripts.txt
  std::unordered_map<std::string, int32_t> txn_tid;
  std::vector<std::vector<int32_t>> tids_per_file; // list of tids as they occur for each file
  std::vector<int32_t> tids;                       // a vector of tids
  int32_t tid = 0;

  std::ofstream ofn(opt.output + "/transcripts.txt"); // master transcripts.txt

  // iterate through each file and populate txn_tid, tids perfile
  for (int i = 0; i < nf; i++)
  {
    tids.clear();
    std::ifstream ifn(opt.files[i] + "/transcripts.txt");
    std::string txn;

    while (ifn >> txn) // while still have transcript data
    {
      auto ok = txn_tid.insert({txn, tid}); // insert transcript and new index assoc with it
      if (ok.second)                        // if the insertion successful
      {
        tids.push_back(tid); // add the index to tids (this is a list of new index)
        ofn << txn << "\n";  // write to file
        tid += 1;            // increment transcript index
      }
      else
      {
        tids.push_back(ok.first->second); // this only hapens when the transcript appears in more than one busfile
      }
    }
    tids_per_file.push_back(tids); // new index
  }

  ofn.close();
  std::cout << "[info] parsed transcripts.txt" << std::endl;

  // all of the ecs are in the header
  // h.ecs is a vector<vector<ints>>
  // the first positional index is the equivalence eid
  // so h.ecs[eid_1] -> returns a vector of tids wrt local indexing
  BUSHeader oh;
  oh.version = BUSFORMAT_VERSION;
  oh.text = "Merged files";

  oh.bclen = vh[0].bclen;
  oh.umilen = vh[0].umilen;

  std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> ecmapinv; // set{tids} (ec) to eid it came from

  for (int i = 0; i < tid; i++)
  {
    oh.ecs.push_back({i});
    ecmapinv.insert({{i}, i});
  }

  std::vector<std::vector<int32_t>> eids_per_file;
  std::vector<int32_t> eids;
  int32_t eid = ecmapinv.size();

  for (int i = 0; i < nf; i++)
  {
    eids.clear();

    BUSHeader h = vh[i];
    const auto &tids = tids_per_file[i]; // new index of tids for that file of the length of that file

    for (const auto &ecs : h.ecs) // ecs is a set of tids std::vector<int32_t>
    {
      std::vector<int32_t> new_ecs(ecs.size());

      // convert tid to new coordinates
      for (int j = 0; j < ecs.size(); j++)
      {
        new_ecs[j] = tids[ecs[j]];
      }

      // check to see if the set exists in ecmapinv
      std::sort(new_ecs.begin(), new_ecs.end());
      new_ecs.erase(std::unique(new_ecs.begin(), new_ecs.end()), new_ecs.end()); // keep only one of the duplicates
      auto it = ecmapinv.find(new_ecs);                                          // see if new_ecs exists
      if (it != ecmapinv.end())
      {
        eid = it->second; // return the eid that it corresponds to
      }
      else
      {
        eid = ecmapinv.size();     // make new eid
        oh.ecs.push_back(new_ecs); // add the set of tids (new ref)
        ecmapinv.insert({new_ecs, eid});
      }
      eids.push_back(eid);
    }
    eids_per_file.push_back(std::move(eids));
  }
  std::cout << "[info] parsed matrix.ec files" << std::endl;
  // generate ecmap from ecmapinv
  std::vector<std::vector<int32_t>> ecmap(ecmapinv.size());
  for (const auto &ec : ecmapinv)
  {
    ecmap[ec.second] = ec.first; // eid -> ecs (set of tids)
  }

  // Process the busfiles
  std::ofstream outf(opt.output + "/merged.bus");
  writeHeader(outf, oh);
  // writeECs(opt.output + "/raw.ec", oh); // prior to reading bus records

  std::priority_queue<BP, std::vector<BP>, std::function<bool(const BP &a, const BP &b)>> pq(pqcmp);
  BUSData bd1, bd2;

  size_t nr = 0, nw = 0;
  // fill the pq
  int bufsize = 100;
  for (int i = 0; i < nf; i++)
  {
    // if we can read data
    if (bf[i].good())
    {
      bf[i].read((char *)&bd1, sizeof(bd1));
      bf[i].read((char *)&bd2, sizeof(bd2));
      assert(bd1.barcode == bd2.barcode && bd1.UMI == bd2.UMI && bd1.ec == bd2.ec && bd1.flags == bd2.flags);

      int32_t prev_rn = bd1.flags;
      while (prev_rn == bd2.flags && bf[i].good())
      {
        pq.push({{bd1, bd2}, i});
        nr += 2;

        bf[i].read((char *)&bd1, sizeof(bd1));
        bf[i].read((char *)&bd2, sizeof(bd2));

        assert(bd1.barcode == bd2.barcode && bd1.UMI == bd2.UMI && bd1.ec == bd2.ec && bd1.flags == bd2.flags);
      }
      // if the flag changes, we want to move the file pointer back by two busrecords
      bf[i].seekg(-sizeof(bd1), std::ios::cur);
      bf[i].seekg(-sizeof(bd1), std::ios::cur);
    }
  }
  // when condition breaks, bd1 and bd2 not pushed to pq, add it in while loop
  // std::cout << "size of pq: " << pq.size() << std::endl;
  // while (!pq.empty())
  // {
  //   print_bd(pq.top().first.first, oh.bclen, oh.umilen);
  //   print_bd(pq.top().first.second, oh.bclen, oh.umilen);
  //   pq.pop();
  // }
  BPair prev_bp, curr_bp;
  int prev_file, curr_file;

  prev_bp = pq.top().first;
  prev_file = pq.top().second;

  prev_bp.first.ec = eids_per_file[prev_file][prev_bp.first.ec];
  prev_bp.second.ec = eids_per_file[prev_file][prev_bp.second.ec];

  pq.pop();

  std::vector<BUSRange> itv_left, itv_right;
  int loop_num = 0;
  while (!pq.empty())
  {
    // std::cout << "- Loop number: " << loop_num << std::endl;
    // std::cout << "      Read number: " << prev_bp.first.flags << std::endl;

    curr_bp = pq.top().first;
    curr_file = pq.top().second;
    curr_bp.first.ec = eids_per_file[curr_file][curr_bp.first.ec];
    curr_bp.second.ec = eids_per_file[curr_file][curr_bp.second.ec];

    pq.pop();

    if (prev_bp.first.flags == curr_bp.first.flags && prev_bp.first.barcode == curr_bp.first.barcode && prev_bp.first.UMI == curr_bp.first.UMI)
    {

      // std::cout << "From file " << prev_file << std::endl;
      // print_bd(prev_bp.first, oh.bclen, oh.umilen);
      // print_bd(prev_bp.second, oh.bclen, oh.umilen);

      itv_left.push_back({prev_bp.first, {prev_bp.first.pad, prev_bp.second.pad}});
      itv_right.push_back({prev_bp.first, {prev_bp.first.pad, prev_bp.second.pad}});
    }
    else
    {

      // std::cout << "From file " << prev_file << std::endl;
      // print_bd(prev_bp.first, oh.bclen, oh.umilen);
      // print_bd(prev_bp.second, oh.bclen, oh.umilen);

      BUSData BR = prev_bp.first;
      itv_left.push_back({prev_bp.first, {prev_bp.first.pad, prev_bp.second.pad}});
      itv_right.push_back({prev_bp.first, {prev_bp.first.pad, prev_bp.second.pad}});

      // 1. sort by lower bound
      std::sort(itv_left.begin(), itv_left.end(), itv_cmp_left);
      // std::cout << "Original intervals: (L, R)\t"
      //           << "(" << itv_left.size() << ", " << itv_right.size() << ")" << std::endl;
      // for (int i = 0; i < itv_left.size(); i++)
      // {
      //   std::cout << itv_left[i].first.ec << ", " << itv_right[i].first.ec << ":\t"
      //             << "[" << itv_left[i].second.first << ", " << itv_left[i].second.second << ")\t"
      //             << "[" << itv_right[i].second.first << ", " << itv_right[i].second.second << ")" << std::endl;
      // }
      // 2. sort by upper bound
      std::sort(itv_right.begin(), itv_right.end(), itv_cmp_right);
      // 3. find elementary sets

      std::pair<std::vector<std::unordered_set<int32_t>>, std::vector<std::pair<int32_t, int32_t>>> thepair = find_elem(itv_left, itv_right);
      std::vector<std::unordered_set<int32_t>> elem_sets = thepair.first;
      std::vector<std::pair<int32_t, int32_t>> bounds = thepair.second;
      // std::cout << "Elementary intervals" << std::endl;
      // for (int i = 0; i < bounds.size(); i++)
      // {
      //   std::cout << "[" << bounds[i].first << ", " << bounds[i].second << ")" << std::endl;
      // }

      // std::cout << "before intersecting, tids, #elem intvs: " << elem_sets.size() << std::endl;
      // for (int i = 0; i < elem_sets.size(); i++)
      // {
      //   for (auto &s : elem_sets[i])
      //   {
      //     for (auto &t : get_tids(oh, ecmapinv, s))
      //     {
      //       std::cout << t << ", ";
      //     }
      //   }
      //   std::cout << std::endl;
      // }

      // 4. intersect tids in the ecs
      // 4a. add the set of tids for a single elem itv
      // only do this if we have elem sets
      if (elem_sets.size() > 0)
      {
        std::vector<int32_t> tids_per_elem, prev_tids;
        for (auto &e : elem_sets[0])
        {
          const auto &tids = oh.ecs[e];
          tids_per_elem.insert(tids_per_elem.end(), tids.begin(), tids.end());
        }
        std::sort(tids_per_elem.begin(), tids_per_elem.end());
        tids_per_elem.erase(std::unique(tids_per_elem.begin(), tids_per_elem.end()), tids_per_elem.end());

        prev_tids = std::move(tids_per_elem);

        // 4b. add the set of tids for the rest of the elem itv
        // if there is more than one elem itv
        if (elem_sets.size() > 1)
        {
          for (int i = 1; i < elem_sets.size(); i++)
          {
            for (auto &e : elem_sets[i])
            {
              const auto &tids = oh.ecs[e];
              tids_per_elem.insert(tids_per_elem.end(), tids.begin(), tids.end());
            }
            std::sort(tids_per_elem.begin(), tids_per_elem.end());
            tids_per_elem.erase(std::unique(tids_per_elem.begin(), tids_per_elem.end()), tids_per_elem.end());

            // Intersect the vectors
            if (tids_per_elem.size())
            {
              prev_tids = intersect_vecs(prev_tids, tids_per_elem);
              tids_per_elem.clear();
            }
          }
        }
        // std::cout << "After intersecting TIDS" << std::endl;
        // for (auto &t : prev_tids)
        // {
        //   std::cout << t << ", ";
        // }
        // std::cout << std::endl;
        if (prev_tids.size()) // if there are tids left over
        {
          // we find the ec associated with prev_tids
          auto findec = ecmapinv.find(prev_tids);
          // if the ec doesnt exist, make a new one
          if (findec == ecmapinv.end())
          {
            BR.ec = ecmapinv.size();
            oh.ecs.push_back(prev_tids);
            ecmapinv.insert({prev_tids, BR.ec});
          }
          else // if it does exist, assign it
          {
            BR.ec = findec->second;
          }
          // write the busrecord
          // std::cout << "ec: " << BR.ec << std::endl;
          BR.count = 1;
          BR.pad = 0;
          outf.write((char *)&BR, sizeof(BR));
          ++nw;
          prev_tids.clear();
        }
        // std::cout << "EC: " << BR.ec << std::endl;
        // std::cout << "##########################" << std::endl;
        tids_per_elem.clear();
      }

      // prev_tids is now a set of tids

      elem_sets.clear();
      bounds.clear();
      itv_left.clear();
      itv_right.clear();
    }
    prev_bp = std::move(curr_bp); // update the bus pair
    prev_file = curr_file;

    // read in more data!
    if (bf[curr_file].good())
    {

      bf[curr_file].read((char *)&bd1, sizeof(bd1));
      bf[curr_file].read((char *)&bd2, sizeof(bd2));
      assert(bd1.barcode == bd2.barcode && bd1.UMI == bd2.UMI && bd1.ec == bd2.ec && bd1.flags == bd2.flags);

      int32_t prev_nr = bd2.flags;
      while (prev_nr == bd2.flags && bf[curr_file].good())
      {
        pq.push({{bd1, bd2}, curr_file});
        nr += 2;

        bf[curr_file].read((char *)&bd1, sizeof(bd1));
        bf[curr_file].read((char *)&bd2, sizeof(bd2));
        assert(bd1.barcode == bd2.barcode && bd1.UMI == bd2.UMI && bd1.ec == bd2.ec && bd1.flags == bd2.flags);
      }
      bf[curr_file].seekg(-sizeof(bd1), std::ios::cur);
      bf[curr_file].seekg(-sizeof(bd2), std::ios::cur);
    }
    loop_num += 1;
  }
  // delete after
  // std::cout << "From file " << prev_file << std::endl;
  // print_bd(prev_bp.first, oh.bclen, oh.umilen);
  // print_bd(prev_bp.second, oh.bclen, oh.umilen);

  BUSData BR = prev_bp.first;
  itv_left.push_back({prev_bp.first, {prev_bp.first.pad, prev_bp.second.pad}});
  itv_right.push_back({prev_bp.first, {prev_bp.first.pad, prev_bp.second.pad}});

  // 1. sort by lower bound
  std::sort(itv_left.begin(), itv_left.end(), itv_cmp_left);
  // std::cout << "Original intervals: (L, R)\t"
  //           << "(" << itv_left.size() << ", " << itv_right.size() << ")" << std::endl;
  // for (int i = 0; i < itv_left.size(); i++)
  // {
  //   std::cout << itv_left[i].first.ec << ", " << itv_right[i].first.ec << ":\t"
  //             << "[" << itv_left[i].second.first << ", " << itv_left[i].second.second << ")\t"
  //             << "[" << itv_right[i].second.first << ", " << itv_right[i].second.second << ")" << std::endl;
  // }
  // 2. sort by upper bound
  std::sort(itv_right.begin(), itv_right.end(), itv_cmp_right);
  // 3. find elementary sets

  std::pair<std::vector<std::unordered_set<int32_t>>, std::vector<std::pair<int32_t, int32_t>>> thepair = find_elem(itv_left, itv_right);
  std::vector<std::unordered_set<int32_t>> elem_sets = thepair.first;
  std::vector<std::pair<int32_t, int32_t>> bounds = thepair.second;
  // std::cout << "Elementary intervals" << std::endl;
  // for (int i = 0; i < bounds.size(); i++)
  // {
  //   std::cout << "[" << bounds[i].first << ", " << bounds[i].second << ")" << std::endl;
  // }

  // std::cout << "before intersecting, tids, #elem intvs: " << elem_sets.size() << std::endl;
  // for (int i = 0; i < elem_sets.size(); i++)
  // {
  //   for (auto &s : elem_sets[i])
  //   {
  //     for (auto &t : get_tids(oh, ecmapinv, s))
  //     {
  //       std::cout << t << ", ";
  //     }
  //   }
  //   std::cout << std::endl;
  // }

  // 4. intersect tids in the ecs
  // 4a. add the set of tids for a single elem itv
  // only do this if we have elem sets
  if (elem_sets.size() > 0)
  {
    std::vector<int32_t> tids_per_elem, prev_tids;
    for (auto &e : elem_sets[0])
    {
      const auto &tids = oh.ecs[e];
      tids_per_elem.insert(tids_per_elem.end(), tids.begin(), tids.end());
    }
    std::sort(tids_per_elem.begin(), tids_per_elem.end());
    tids_per_elem.erase(std::unique(tids_per_elem.begin(), tids_per_elem.end()), tids_per_elem.end());

    prev_tids = std::move(tids_per_elem);

    // 4b. add the set of tids for the rest of the elem itv
    // if there is more than one elem itv
    if (elem_sets.size() > 1)
    {
      for (int i = 1; i < elem_sets.size(); i++)
      {
        for (auto &e : elem_sets[i])
        {
          const auto &tids = oh.ecs[e];
          tids_per_elem.insert(tids_per_elem.end(), tids.begin(), tids.end());
        }
        std::sort(tids_per_elem.begin(), tids_per_elem.end());
        tids_per_elem.erase(std::unique(tids_per_elem.begin(), tids_per_elem.end()), tids_per_elem.end());

        // Intersect the vectors
        if (tids_per_elem.size())
        {
          prev_tids = intersect_vecs(prev_tids, tids_per_elem);
          tids_per_elem.clear();
        }
      }
    }
    // std::cout << "After intersecting TIDS" << std::endl;
    // for (auto &t : prev_tids)
    // {
    //   std::cout << t << ", ";
    // }
    // std::cout << std::endl;
    if (prev_tids.size()) // if there are tids left over
    {
      // we find the ec associated with prev_tids
      auto findec = ecmapinv.find(prev_tids);
      // if the ec doesnt exist, make a new one
      if (findec == ecmapinv.end())
      {
        BR.ec = ecmapinv.size();
        oh.ecs.push_back(prev_tids);
        ecmapinv.insert({prev_tids, BR.ec});
      }
      else // if it does exist, assign it
      {
        BR.ec = findec->second;
      }
      // write the busrecord
      // std::cout << "ec: " << BR.ec << std::endl;
      BR.count = 1;
      BR.pad = 0;
      outf.write((char *)&BR, sizeof(BR));
      ++nw;
      prev_tids.clear();
    }
    // std::cout << "EC: " << BR.ec << std::endl;
    // std::cout << "##########################" << std::endl;
  }

  // prev_tids is now a set of tids

  elem_sets.clear();
  itv_left.clear();
  itv_right.clear();
  // close all of the busfiles
  for (int i = 0; i < nf; i++)
  {
    bf[i].close();
  }
  // write the matrix.ec
  writeECs(opt.output + "/matrix.ec", oh);
  std::cerr << "bus records read:    " << nr << std::endl;
  std::cerr << "bus records written: " << nw << std::endl;
}
