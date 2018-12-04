#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>
#include <thread>
#include <atomic>
#include <unordered_map>
#include <functional>

#include "Common.hpp"
#include "BUSData.h"

int my_mkdir(const char *path, mode_t mode) {
  #ifdef _WIN64
  return mkdir(path);
  #else
  return mkdir(path,mode);
  #endif
}

bool checkFileExists(const std::string &fn) {
  struct stat stFileInfo;
  auto intStat = stat(fn.c_str(), &stFileInfo);
  return intStat == 0;
}

bool checkDirectoryExists(const std::string &fn) {
  struct stat stFileInfo;
  auto intStat = stat(fn.c_str(), &stFileInfo);
  return intStat == 0 && S_ISDIR(stFileInfo.st_mode);
}


struct SortedVectorHasher {
  size_t operator()(const std::vector<int32_t>& v) const {
    uint64_t r = 0;
    int i=0;
    for (auto x : v) {
      uint64_t t = std::hash<int32_t>{}(x);
      t = (x>>i) | (x<<(64-i));
      r = r ^ t;
      i = (i+1)%64;
    }
    return r;
  }
};



void parse_ProgramOptions_sort(int argc, char **argv, Bustools_opt& opt) {

  const char* opt_string = "t:o:";

  static struct option long_options[] = {

    {"threads",         required_argument,  0, 't'},
    {"output",          required_argument,  0, 'o'},
    {0,                 0,                  0,  0 }
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

    switch (c) {

    case 't':
      opt.threads = atoi(optarg);
      break;    
    case 'o':
      opt.output = optarg;
      break;
    default:
      break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  while (optind < argc) opt.files.push_back(argv[optind++]);
}

void parse_ProgramOptions_merge(int argc, char **argv, Bustools_opt& opt) {
   const char* opt_string = "o:";

  static struct option long_options[] = {
    {"output",          required_argument,  0, 'o'},
    {0,                 0,                  0,  0 }
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

    switch (c) {
    case 'o':
      opt.output = optarg;
      break;
    default:
      break;
    }
  }

  while (optind < argc) opt.files.push_back(argv[optind++]);
}

void parse_ProgramOptions_dump(int argc, char **argv, Bustools_opt& opt) {

  const char* opt_string = "o:";

  static struct option long_options[] = {
    {"output",          required_argument,  0, 'o'},
    {0,                 0,                  0,  0 }
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

    switch (c) {
    case 'o':
      opt.output = optarg;
      break;
    default:
      break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  while (optind < argc) opt.files.push_back(argv[optind++]);
}



bool check_ProgramOptions_sort(Bustools_opt& opt) {

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

  if (opt.output.empty()) {
    std::cerr << "Error missing output file" << std::endl;
    ret = false;
  } 


  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  } else {
    for (const auto& it : opt.files) {  
      if (!checkFileExists(it)) {
        std::cerr << "Error: File not found, " << it << std::endl;
        ret = false;
      }
    }
  }

  return ret;
}


bool check_ProgramOptions_merge(Bustools_opt& opt) {
  bool ret = true;
  
  if (opt.output.empty()) {
    std::cerr << "Error missing output directory" << std::endl;
    ret = false;
  } 

  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing BUS input directories" << std::endl;
    ret = false;
  } else {
    for (const auto& it : opt.files) {
      if (!checkDirectoryExists(it)) {
        std::cerr << "Error: directory " << it << " does not exist or is not a directory" << std::endl;
        ret = false;
      }
    }
    if (ret) {
      // check for output.bus and matrix.ec in each of the directories
      for (const auto &it : opt.files) {
        if (!checkFileExists(it + "/output.bus")) {
          std::cerr << "Error: file " << it << "/output.bus not found" << std::endl;
        }

        if (!checkFileExists(it + "/matrix.ec")) {
          std::cerr << "Error: file " << it << "/matrix.ec not found" << std::endl;
        }
      }
    }

    // check if output directory exists or if we can create it
    struct stat stFileInfo;
    auto intStat = stat(opt.output.c_str(), &stFileInfo);
    if (intStat == 0) {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode)) {
        std::cerr << "Error: file " << opt.output << " exists and is not a directory" << std::endl;
        ret = false;
      } 
    } else {
      // create directory
      if (my_mkdir(opt.output.c_str(), 0777) == -1) {
        std::cerr << "Error: could not create directory " << opt.output << std::endl;
        ret = false;
      }
    }
  }

  return ret;

}

bool check_ProgramOptions_dump(Bustools_opt& opt) {
  bool ret = true;

  if (opt.output.empty()) {
    std::cerr << "Error missing output file" << std::endl;
    ret = false;
  } 


  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  } else {
    for (const auto& it : opt.files) {  
      if (!checkFileExists(it)) {
        std::cerr << "Error: File not found, " << it << std::endl;
        ret = false;
      }
    }
  }

  return ret;
}



void Bustools_Usage() {
  std::cout << "Usage: bustools <CMD> [arguments] .." << std::endl << std::endl
  << "Where <CMD> can be one of: " << std::endl << std::endl
  << "sort            Sort bus file by barcodes and UMI" << std::endl
  << "text            Output as tab separated text file" << std::endl << std::endl
  << "Running bustools <CMD> without arguments prints usage information for <CMD>"
  << std::endl << std::endl;
}



void Bustools_sort_Usage() {
  std::cout << "Usage: bustools sort [options] bus-files" << std::endl << std::endl
  << "Options: " << std::endl
  << "-t, --threads         Number of threads to use" << std::endl
  << "-o, --output          File for sorted output" << std::endl
  << std::endl;
}

void Bustools_merge_Usage() {
  std::cout << "Usage: bustools merge [options] directories" << std::endl << std::endl
  << "Options: " << std::endl
  << "-t, --threads         Number of threads to use" << std::endl
  << "-o, --output          Directory for merged output" << std::endl
  << std::endl;
}

void Bustools_dump_Usage() {
  std::cout << "Usage: bustools text [options] bus-files" << std::endl << std::endl
  << "Options: " << std::endl
  << "-o, --output          File for text output" << std::endl
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






int main(int argc, char **argv) {

  if (argc < 2) {
    // Print error message, function?
    Bustools_Usage();
    exit(1);
  } else {
    bool disp_help = argc == 2;
    std::string cmd(argv[1]);
    Bustools_opt opt;


    if (cmd == "sort") {
      if (disp_help) {
        Bustools_sort_Usage();
        exit(0);        
      }
      parse_ProgramOptions_sort(argc-1, argv+1, opt);
      if (check_ProgramOptions_sort(opt)) { //Program options are valid
        BUSHeader h;
        std::vector<BUSData> b;
        size_t N = 100000;
        BUSData* p = new BUSData[N];
        char magic[4];
        uint32_t version = 0;
        for (const auto& infn : opt.files) { 
          std::ifstream inf(infn.c_str(), std::ios::binary);
          parseHeader(inf, h);

          int rc = 1;
          while (true) {
            inf.read((char*)p, N*sizeof(BUSData));
            size_t rc = inf.gcount() / sizeof(BUSData);
            if (rc == 0) {
              break;
            }
            b.insert(b.end(), p, p+rc);
          }
        }
        delete[] p; p = nullptr;
        std::cerr << "Read in " << b.size() << " number of busrecords" << std::endl;
        std::sort(b.begin(), b.end(), [&](const BUSData& a, const BUSData &b) 
                                        {
                                          if (a.barcode == b.barcode) {
                                          if (a.UMI == b.UMI) {
                                            return a.ec < b.ec;
                                          } else {
                                           return a.UMI < b.UMI;
                                          } 
                                         } else { 
                                           return a.barcode < b.barcode;
                                         }});
        std::cerr << "All sorted" << std::endl;

        std::ofstream busf_out;
        busf_out.open(opt.output , std::ios::out | std::ios::binary);
        writeHeader(busf_out, h);

        size_t n = b.size();
        for (size_t i = 0; i < n; ) {
          size_t j = i+1;
          uint32_t c = b[i].count;
          auto ec = b[i].ec;          
          for (; j < n; j++) {
            if (b[i].barcode != b[j].barcode || b[i].UMI != b[j].UMI || b[i].ec != b[j].ec) {
              break;
            }
            c += b[j].count;
          }
          // merge identical things
          b[i].count = c;
          busf_out.write((char*)(&(b[i])), sizeof(b[i]));
          // increment
          i = j;
        }
        busf_out.close();    
      } else {
        Bustools_sort_Usage();
        exit(1);
      }
    } else if (cmd == "merge") {
      if (disp_help) {
        Bustools_merge_Usage();
        exit(0);
      }
      parse_ProgramOptions_merge(argc-1, argv+1, opt);
      if (check_ProgramOptions_merge(opt)) {
        // first parse all headers
        std::vector<BUSHeader> vh;
        // TODO: check for compatible headers, version numbers umi and bclen

        for (const auto& infn : opt.files) {
          std::ifstream inf((infn + "/output.bus").c_str(), std::ios::binary);
          BUSHeader h;
          parseHeader(inf, h);
          inf.close();
          
          parseECs(infn + "/matrix.ec", h);
          vh.push_back(std::move(h));
        }

        // create master ec
        BUSHeader oh;
        oh.version = BUSFORMAT_VERSION;
        oh.text = "Merged files from BUStools";
        //TODO: parse the transcripts file, check that they are identical and merge.
        oh.bclen = vh[0].bclen;
        oh.umilen = vh[0].umilen;
        std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> ecmapinv;
        std::vector<std::vector<int32_t>> ectrans;        
        std::vector<int32_t> ctrans;
        
        oh.ecs = vh[0].ecs; // copy operator

        for (int32_t ec = 0; ec < oh.ecs.size(); ec++) {
          ctrans.push_back(ec);
          const auto &v = oh.ecs[ec];
          ecmapinv.insert({v, ec});
        }
        ectrans.push_back(std::move(ctrans));
        
        for (int i = 1; i < opt.files.size(); i++) {
          ctrans.clear();
          // merge the rest of the ecs
          int j = -1;
          for (const auto &v : vh[i].ecs) {
            j++;
            int32_t ec = -1;
            auto it = ecmapinv.find(v);
            if (it != ecmapinv.end()) {
              ec = it->second;              
            } else {
              ec = ecmapinv.size();
              oh.ecs.push_back(v); // copy
              ecmapinv.insert({v,ec});
            }
            ctrans.push_back(ec);
          }
          ectrans.push_back(ctrans);
        }

        // now create a single output file
        writeECs(opt.output + "/matrix.ec", oh);
        std::ofstream outf(opt.output + "/output.bus");
        writeHeader(outf, oh);


        size_t N = 100000;
        BUSData* p = new BUSData[N];
        size_t nr = 0;
        for (int i = 0; i < opt.files.size(); i++) {
          // open busfile and parse header
          BUSHeader h;
          const auto &ctrans = ectrans[i];
          std::ifstream inf((opt.files[i] + "/output.bus").c_str(), std::ios::binary);
          parseHeader(inf, h);
          // now read all records and translate the ecs
          while (true) {
            inf.read((char*)p, N*sizeof(BUSData));
            size_t rc = inf.gcount() / sizeof(BUSData);
            if (rc == 0) {
              break;
            }
            nr += rc;
            for (size_t i = 0; i < rc; i++) {
              auto &b = p[i];
              b.ec = ctrans[b.ec]; // modify the ec              
            }
            outf.write((char*)p, rc*sizeof(BUSData));
          }
          inf.close();
        }
        outf.close();
      } else {
        Bustools_merge_Usage();
        exit(1);
      }

    } else if (cmd == "dump" || cmd == "text") {
      if (disp_help) {
        Bustools_dump_Usage();
        exit(0);        
      }
      parse_ProgramOptions_dump(argc-1, argv+1, opt);
      if (check_ProgramOptions_dump(opt)) { //Program options are valid
        BUSHeader h;
        size_t nr = 0;
        size_t N = 100000;
        BUSData* p = new BUSData[N];
        std::ofstream of;
        of.open(opt.output); 
        char magic[4];      
        uint32_t version = 0;
        for (const auto& infn : opt.files) { 
          std::ifstream inf(infn.c_str(), std::ios::binary);
          parseHeader(inf, h);
          uint32_t bclen = h.bclen;
          uint32_t umilen = h.umilen;
          int rc = 0;
          while (true) {
            inf.read((char*)p, N*sizeof(BUSData));
            size_t rc = inf.gcount() / sizeof(BUSData);
            if (rc == 0) {
              break;
            }
            nr += rc;
            for (size_t i = 0; i < rc; i++) {
              of << binaryToString(p[i].barcode, bclen) << "\t" << binaryToString(p[i].UMI,umilen) << "\t" << p[i].ec << "\t" << p[i].count << "\n";        
            }
          }
        }
        delete[] p; p = nullptr;
        of.close();
        std::cerr << "Read in " << nr << " number of busrecords" << std::endl;
      } else {
        Bustools_dump_Usage();
        exit(1);
      }
    } else {
      std::cerr << "Error: invalid command " << cmd << std::endl;
      Bustools_Usage();      
    }

  }
}
