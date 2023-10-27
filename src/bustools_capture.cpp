#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_capture.h"

void bustools_capture(Bustools_opt &opt) {
  BUSHeader h;

  std::unordered_set<uint64_t> captures;
  std::vector<std::vector<int32_t>> ecmap;
  u_map_<std::vector<int32_t>, int32_t, SortedVectorHasher> ecmapinv;
  bool capture_prefixes = false;

  if (opt.type == CAPTURE_TX) {
    // parse ecmap and capture list
    u_map_<std::string, int32_t> txnames;
    std::cerr << "Parsing transcripts .. "; std::cerr.flush();
    parseTranscripts(opt.count_txp, txnames);
    std::cerr << "done" << std::endl;
    std::cerr << "Parsing ECs .. "; std::cerr.flush();
    parseECs(opt.count_ecs, h);
    std::cerr << "done" << std::endl;
    ecmap = h.ecs; // copy

    ecmapinv.reserve(ecmap.size());
    for (int32_t ec = 0; ec < ecmap.size(); ec++) {
      ecmapinv.insert({ecmap[ec], ec});
    }

    std::cerr << "Parsing capture list .. "; std::cerr.flush();
    parseTxCaptureList(opt.capture, txnames, captures);
    std::cerr << "done" << std::endl;
  } else if (opt.type == CAPTURE_UMI || opt.type == CAPTURE_BC) {
    parseBcUmiCaptureList(opt.capture, captures);
    if (opt.type == CAPTURE_BC && captures.count(std::numeric_limits<uint64_t>::max()) > 0) {
      capture_prefixes = true;
    }
  } else if (opt.type == CAPTURE_F) {
    parseFlagsCaptureList(opt.capture, captures);
  } else { // Should never happen
    std::cerr << "ERROR: Unknown capture type" << std::endl;
    exit(1);
  }

  bool outheader_written = false;

  std::string output = opt.output;
  if (opt.filter) {
    output += ".bus";
  }

  std::streambuf *buf = nullptr;
  std::ofstream of;
  if (!opt.stream_out) {
    of.open(output); 
    buf = of.rdbuf();
  } else {
    buf = std::cout.rdbuf();
  }
  std::ostream o(buf);

  size_t nr = 0, nw = 0;
  size_t N = 100000;
  BUSData* p = new BUSData[N];
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
    if (h.bclen == 32) capture_prefixes = false;
    uint64_t len_mask = ((1ULL << (2*h.bclen)) - 1);

    if (!outheader_written) {
      writeHeader(o, h);
      outheader_written = true;
    }

    while(true) {
      in.read((char*)p, N*sizeof(BUSData));
      size_t rc = in.gcount() / sizeof(BUSData);
      if (rc == 0) {
        break;
      }
      nr += rc;

      for (size_t i = 0; i < rc; i++) {
        bd = p[i];
        bool capt = false;

        if (opt.type == CAPTURE_TX) {
          if (bd.ec < 0 || bd.ec > ecmap.size()) {
            continue;
          }
          for (auto x : ecmap[bd.ec]) {
            if (captures.count((uint64_t) x) > 0) {
              capt = true;
              break;
            }
          }

        } else if (opt.type == CAPTURE_BC) {
          if (capture_prefixes) {
            uint64_t bitmask = (1ULL << (2*(32-bd.bclen))) - 1;
            uint64_t potential_prefix_barcode = (bd.barcode >> (2*bd.bclen)) & bitmask;
            // Now we have the sequence preceding the actual cell barcode, let's find it in the captures set
            capt = captures.count((static_cast<uint64_t>(1) << 63) | potential_prefix_barcode) > 0; // If prefix exists in the "capture prefix" set
            if (!capt) capt = captures.count(bd.barcode & len_mask) > 0; // If not, then check the barcode as-is
          } else {
            capt = captures.count(bd.barcode & len_mask) > 0;
          }
        } else if (opt.type == CAPTURE_UMI) {
          capt = captures.count(bd.UMI) > 0;
        } else if (opt.type == CAPTURE_F) {
          capt = captures.count(bd.flags) > 0;
        } else { // Should never happen
          std::cerr << "error: unknown capture type" << std::endl;
          exit(1);
        }
        
        if (capt != opt.complement) {
          if (opt.filter) { // modify the ec
            std::vector<int32_t> v;
            for (auto x : ecmap[bd.ec]) {
              if (captures.count((uint64_t) x) > 0) {
                v.push_back(x);
              }
            }
            
            if (v.empty()) {
              continue; // should never happen
            } else {
              std::sort(v.begin(), v.end());                  
            }

            auto it = ecmapinv.find(v);
            if (it == ecmapinv.end()) {
              // create new ec;
              int32_t ec = ecmap.size();
              ecmap.push_back(v);
              ecmapinv.insert({v,ec});
              bd.ec = ec;
            } else {
              bd.ec = it->second;
            }
          }

          ++nw;
          o.write((char *) &bd, sizeof(bd));
        }
      }
    }
    if (!opt.stream_in) {
      inf.close();
    }
  }

  if (opt.filter) {
    h.ecs = ecmap; // modified map
    // TODO: trim down the ecs for the capture list
    writeECs(opt.output + ".ec", h);
  }

  if (of.is_open()) {
    of.close();
  }

  std::cerr << "Read in " << nr << " BUS records, wrote " << nw << " BUS records" << std::endl;
}
