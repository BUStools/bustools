#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>
#include <thread>
#include <atomic>

#include "Common.hpp"
#include "BUSData.h"

void parse_ProgramOptions_barcode(int argc, char **argv, Bustools_opt& opt) {

  const char* opt_string = "t:";

  static struct option long_options[] = {

    {"threads",         required_argument,  0, 't'},
    {0,                 0,                  0,  0 }
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

    switch (c) {

    case 't':
      opt.threads = atoi(optarg);
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
  }
  else if (opt.threads > max_threads) {

    std::cerr << "Warning: Number of threads cannot be greater than or equal to " << max_threads 
    << ". Setting number of threads to " << max_threads << std::endl;
    opt.threads = max_threads;
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
  << std::endl;
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
        for (size_t i = 0; i < n; i++) {
          size_t j = i+1;
          for (; j < n; j++) {
            if (b[i].barcode != b[j].barcode || b[i].UMI != b[j].UMI) {
              break;
            }
          }
          // b[i,j) have same barcode and UMI and is maximal
          size_t dz = j - i;
          if (dz > 1) {
            std::cout << binaryToString(b[i].barcode, 16) << "\t" << binaryToString(b[i].UMI,10);
            for (size_t k = i; k < j; k++) {
              std::cout << "\t" << b[k].ec;
            }
            std::cout << "\n";
          }
        }

      } else {
        Bustools_Usage();
        exit(1);
      }
    }
  }
}
