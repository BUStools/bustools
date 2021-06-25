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

inline std::vector<int32_t> get_tids(const BUSHeader &oh, const std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> &ecmapinv, const int32_t &eid)
{

    std::vector<int32_t> tids = oh.ecs[eid];

    std::sort(tids.begin(), tids.end());
    tids.erase(std::unique(tids.begin(), tids.end()), tids.end());

    return tids;
}

inline void print_bd(const BUSData &bd, const size_t bclen, const size_t umilen)
{
    std::cerr << binaryToString(bd.barcode, bclen) << "\t" << binaryToString(bd.UMI, umilen) << "\t" << bd.ec << "\t" << bd.count << "\t" << bd.flags << "\t" << bd.pad << std::endl;
}

// vocab
// txn: string representing transcript str
// tid: number indexing the transcript in the file
// ecs: the set of tids (set of transcript indexes)
// eid: the number indexing the equivalence class
// keep in mind this could be wrt original file or
// merged file

void bustools_mash(const Bustools_opt &opt)
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
    std::cerr << "[info] parsed output.bus files" << std::endl;

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
    std::cerr << "[info] parsed transcripts.txt" << std::endl;

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

    for (int32_t i = 0; i < tid; i++)
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
            for (int32_t j = 0; j < ecs.size(); j++)
            {
                new_ecs[j] = tids[ecs[j]];
            }

            // check to see if the set exists in ecmapinv
            if (new_ecs.size() == 1)
            {
                eid = new_ecs[0];
                if (eid > ecmapinv.size())
                {
                    std::cerr << "[warn] transcript not found" << std::endl;
                }
            }
            else
            {
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
            }

            eids.push_back(eid);
        }
        eids_per_file.push_back(std::move(eids));
    }
    std::cerr << "[info] parsed matrix.ec files" << std::endl;
    // generate ecmap from ecmapinv
    // std::vector<std::vector<int32_t>> ecmap(ecmapinv.size());
    // for (const auto &ec : ecmapinv)
    // {
    //     ecmap[ec.second] = ec.first; // eid -> ecs (set of tids)
    // }

    // Process the busfiles
    std::ofstream outf(opt.output + "/mashed.bus");
    writeHeader(outf, oh);
    // writeECs(opt.output + "/raw.ec", oh); // prior to reading bus records

    BUSData bd, new_bd;
    int32_t nr = 0, nw = 0;

    int32_t N = 1024;
    std::queue<BUSData> queue[nf];

    for (int i = 0; i < nf; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (bf[i].good())
            {
                bf[i].read((char *)&bd, sizeof(bd));
                queue[i].push(bd);
                nr++;
            }
        }
    }

    for (int i = 0; i < nf; i++)
    {
        // if we can read a bus record
        while (true)
        {
            if (queue[i].empty())
            {
                break;
            }

            new_bd = queue[i].front();
            queue[i].pop();

            // print_bd(new_bd, oh.bclen, oh.umilen);
            new_bd.ec = eids_per_file[i][new_bd.ec];
            // print_bd(new_bd, oh.bclen, oh.umilen);
            // std::cout << "------------------------" << std::endl;
            outf.write((char *)&new_bd, sizeof(new_bd));
            nw++;

            if (bf[i].good())
            {
                bf[i].read((char *)&bd, sizeof(bd));
                nr++;
                queue[i].push(bd);
            }
        }
    }

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