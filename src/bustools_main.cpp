#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <atomic>

#include "Common.hpp"
#include "BUSData.h"

void parse_ProgramOptions_barcode(int argc, char **argv, Bustools_opt& opt) {

  const char* opt_string = "t:e:";

  static struct option long_options[] = {

    {"threads",         required_argument,  0, 't'},
    {"ecs",             required_argument,  0, 'e'},
    {0,                 0,                  0,  0 }
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

    switch (c) {

    case 't':
      opt.threads = atoi(optarg);
      break;    
    case 'e':
      opt.ecf = optarg;
      break;
    default:
      break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  while (optind < argc) opt.files.push_back(argv[optind++]);
}

bool check_ProgramOptions_barcode(Bustools_opt& opt) {

  bool ret = true;

  size_t max_threads = std::thread::hardware_concurrency();

  if (opt.threads <= 0) {

    std::cerr << "Error: Number of threads cannot be less than or equal to 0" << std::endl;
    ret = false;
  } else if (opt.threads > max_threads) {

    std::cerr << "Warning: Number of threads cannot be greater than or equal to " << max_threads 
    << ". Setting number of threads to " << max_threads << std::endl;
    opt.threads = max_threads;
  }

  if (opt.ecf.empty()) {
    std::cerr << "Error missing EC file" << std::endl;
  } else {
    struct stat stFileInfo;
    int intStat;    
    intStat = stat(opt.ecf.c_str(), &stFileInfo);
    if (intStat != 0) {
      std::cerr << "Error: EC file not found " << opt.ecf << std::endl;
      ret = false;
    }
  }


  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    int intStat;

    for (const auto& it : opt.files) {  
      intStat = stat(it.c_str(), &stFileInfo);

      if (intStat != 0) {
        std::cerr << "Error: File not found, " << it << std::endl;
        ret = false;
      }
    }
  }

  return ret;
}



void Bustools_Usage() {
  std::cout << "Usage: bustools barcodes [options] fasta-files" << std::endl << std::endl
  << "Options: " << std::endl
  << "-t, --threads         Number of threads to use" << std::endl
  << "-e, --ecf             EC index file" << std::endl
  << std::endl;
}


std::vector<int32_t> intersect(std::vector<int32_t> &u, std::vector<int32_t> &v) {
  std::vector<int32_t> res;
  auto a = u.begin();
  auto b = v.begin();

  while (a != u.end() && b != v.end()) {
    if (*a < *b) {
      ++a;
    } else if (*b < *a) {
      ++b;
    } else {
      // match
      res.push_back(*a);
      ++a;
      ++b;
    }
  }
  return std::move(res);
}
std::vector<std::vector<int32_t>> readEC(const std::string &fn) {
  std::vector<std::vector<int32_t>> ecs;  
  std::ifstream inf(fn.c_str());
  std::string line, t;
  line.reserve(10000);
  

  std::vector<int32_t> c;
  
  int i = 0;
  while (std::getline(inf, line)) {       
    c.clear();
    int ec = -1;
    if (line.size() == 0) {
      continue;
    }
    std::stringstream ss(line);
    ss >> ec;
    assert(ec == i);
    while (std::getline(ss, t, ',')) {
      c.push_back(std::stoi(t));
    }

    ecs.push_back(std::move(c));
    i++;
  }
  return ecs;
}


int main(int argc, char **argv) {

  if (argc < 2) {
    // Print error message, function?
    Bustools_Usage();
    exit(1);
  } else {

    std::string cmd(argv[1]);
    Bustools_opt opt;


    if (cmd == "barcodes") {
      parse_ProgramOptions_barcode(argc-1, argv+1, opt);
      if (check_ProgramOptions_barcode(opt)) { //Program options are valid
        //TODO refactor into separate method

        auto ecs = readEC(opt.ecf);

        std::vector<BUSData> b;
        size_t N = 100000;
        BUSData* p = new BUSData[N];
        for (const auto& infn : opt.files) { 
          std::ifstream inf(infn.c_str(), std::ios::binary);
          int rc = 1;
          while (true) {
            inf.read((char*)p, N*sizeof(BUSData));
            size_t rc = inf.gcount() / sizeof(BUSData);
            if (rc == 0) {
              break;
            }
            std::cout << rc << std::endl;
            b.insert(b.end(), p, p+rc);
          }
        }
        delete[] p; p = nullptr;
        std::cout << "Read in " << b.size() << " number of busrecords" << std::endl;
        std::sort(b.begin(), b.end(), [&](const BUSData& a, const BUSData &b) 
                                        {if (a.barcode == b.barcode) {
                                           return a.UMI < b.UMI; }
                                         else { 
                                           return a.barcode < b.barcode;
                                         }});
        std::cout << "All sorted" << std::endl;

        size_t n = b.size();
        std::vector<int32_t> tmp;
        tmp.resize(1000);
        for (size_t i = 0; i < n; ) {
          size_t j = i+1;
          for (; j < n; j++) {
            if (b[i].barcode != b[j].barcode || b[i].UMI != b[j].UMI) {
              break;
            }
          }
          // b[i,j) have same barcode and UMI and is maximal
          size_t dz = j - i;          
          if (dz > 1) {
            tmp.clear();
            for (size_t k = i; k < j; k++) {
              tmp.push_back(b[k].ec);
            }
            std::sort(tmp.begin(), tmp.end());
            int32_t ec = tmp[0];
            std::vector<int32_t> v = ecs[ec];            
            for (size_t k = 1; k < tmp.size(); k++) {
              if (tmp[k] != ec) {
                ec = tmp[k];
                v = intersect(v, ecs[ec]);
              }
              if (v.empty()) {
                break;
              }
            }
            if (v.empty()) {
              std::cout << binaryToString(b[i].barcode, 16) << "\t" << binaryToString(b[i].UMI,10);
              for (size_t k = i; k < j; k++) {
                std::cout << "\t" << b[k].ec;
              }
              std::cout << "\n";
            }
          }
          i = j;
        }


        // todo read 
      } else {
        Bustools_Usage();
        exit(1);
      }
    }
  }
}
