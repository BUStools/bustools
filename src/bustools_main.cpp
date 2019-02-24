#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>
#include <thread>
#include <atomic>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <sstream>

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

  const char* opt_string = "t:o:p";

  static struct option long_options[] = {
    {"threads",         required_argument,  0, 't'},
    {"output",          required_argument,  0, 'o'},
    {"pipe",            no_argument, 0, 'p'},
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
    case 'p':
      opt.stream_out = true;
      break;
    default:
      break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  while (optind < argc) opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-") {
    opt.stream_in = true;
  }
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

void parse_ProgramOptions_count(int argc, char **argv, Bustools_opt& opt) {
  const char* opt_string = "o:g:e:t:";
  int gene_flag = 0;
  static struct option long_options[] = {
    {"output",          required_argument,  0, 'o'},
    {"genemap",          required_argument,  0, 'g'},
    {"ecmap",          required_argument,  0, 'e'},
    {"txnames",          required_argument,  0, 't'},
    {"genecounts", no_argument, &gene_flag, 1},
    {0,                 0,                  0,  0 }
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

    switch (c) {
    case 'o':
      opt.output = optarg;
      break;
    case 'g':
      opt.count_genes = optarg;
      break;
    case 't':
      opt.count_txp = optarg;
      break;
    case 'e':
      opt.count_ecs = optarg;
      break;
    default:
      break;
    }
  }
  if (gene_flag) {
    opt.count_collapse = true;
  }

  while (optind < argc) opt.files.push_back(argv[optind++]);

  if (opt.files.size() == 1 && opt.files[0] == "-") {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_dump(int argc, char **argv, Bustools_opt& opt) {

  const char* opt_string = "o:p";

  static struct option long_options[] = {
    {"output",          required_argument,  0, 'o'},
    {"pipe",            no_argument, 0, 'p'},
    {0,                 0,                  0,  0 }
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

    switch (c) {
    case 'o':
      opt.output = optarg;
      break;
    case 'p':
      opt.stream_out = true;
    default:
      break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  while (optind < argc) opt.files.push_back(argv[optind++]);

  if (opt.files.size() == 1 && opt.files[0] == "-") {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_correct(int argc, char **argv, Bustools_opt& opt) {

  const char* opt_string = "o:w:p";
  static struct option long_options[] = {
    {"output",          required_argument,  0, 'o'},
    {"whitelist",       required_argument,  0, 'w'},
    {"pipe",            no_argument, 0, 'p'},
    {0,                 0,                  0,  0 }
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

    switch (c) {
    case 'o':
      opt.output = optarg;
      break;
    case 'w':
      opt.whitelist = optarg;
      break;
    case 'p':
      opt.stream_out = true;
      break;
    default:
      break;
    }
  }


  //hard code options for now
  opt.ec_d = 1;
  opt.ec_dmin = 3;

  // all other arguments are fast[a/q] files to be read
  while (optind < argc) opt.files.push_back(argv[optind++]);

  if (opt.files.size() == 1 && opt.files[0] == "-") {
    opt.stream_in = true;
  }
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

  if (!opt.stream_out && opt.output.empty()) {
    std::cerr << "Error missing output file" << std::endl;
    ret = false;
  } 


  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  } else if (!opt.stream_in) {
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

  if (!opt.stream_out && opt.output.empty()) {
    std::cerr << "Error missing output file" << std::endl;
    ret = false;
  } 


  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  } else if (!opt.stream_in) {    
    for (const auto& it : opt.files) {  
      if (!checkFileExists(it)) {
        std::cerr << "Error: File not found, " << it << std::endl;
        ret = false;
      }
    }
  }

  return ret;
}

bool check_ProgramOptions_correct(Bustools_opt& opt) {
  bool ret = true;

  if (!opt.stream_out && opt.output.empty()) {
    std::cerr << "Error: Missing output file" << std::endl;
    ret = false;
  } 


  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  } else {
    if (!opt.stream_in) {
      for (const auto& it : opt.files) {  
        if (!checkFileExists(it)) {
          std::cerr << "Error: File not found, " << it << std::endl;
          ret = false;
        }
      }
    }
  }

  if (opt.whitelist.size() == 0) {
    std::cerr << "Error: Missing whitelist file" << std::endl;
    ret = false;
  } else {
    if (!checkFileExists(opt.whitelist)) {
      std::cerr << "Error: File not found " << opt.whitelist << std::endl;
      ret = false;
    }
  }

  return ret;
}

bool check_ProgramOptions_count(Bustools_opt& opt) {
  bool ret = true;

  if (opt.output.empty()) {
    std::cerr << "Error: Missing output file" << std::endl;
    ret = false;
  } 


  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  } else {
    if (!opt.stream_in) {
      for (const auto& it : opt.files) {  
        if (!checkFileExists(it)) {
          std::cerr << "Error: File not found, " << it << std::endl;
          ret = false;
        }
      }
    }
  }

  if (opt.count_genes.size() == 0) {
    std::cerr << "Error: missing gene mapping file" << std::endl;
  } else {
    if (!checkFileExists(opt.count_genes)) {
      std::cerr << "Error: File not found " << opt.count_genes << std::endl;
      ret = false;
    }
  }

  if (opt.count_ecs.size() == 0) {
    std::cerr << "Error: missing equialence class mapping file" << std::endl;
  } else {
    if (!checkFileExists(opt.count_genes)) {
      std::cerr << "Error: File not found " << opt.count_ecs << std::endl;
      ret = false;
    }
  }

  if (opt.count_txp.size() == 0) {
    std::cerr << "Error: missing transcript name file" << std::endl;
  } else {
    if (!checkFileExists(opt.count_genes)) {
      std::cerr << "Error: File not found " << opt.count_txp << std::endl;
      ret = false;
    }
  }

  return ret;
}


void Bustools_Usage() {
  std::cout << "Usage: bustools <CMD> [arguments] .." << std::endl << std::endl
  << "Where <CMD> can be one of: " << std::endl << std::endl
  << "sort            Sort bus file by barcodes and UMI" << std::endl
  << "text            Output as tab separated text file" << std::endl 
  << "merge           Merge bus files from same experiment" << std::endl
  << "correct         Error correct bus files" << std::endl
  << "count           Generate count matrices from bus file" << std::endl
  << std::endl
  << "Running bustools <CMD> without arguments prints usage information for <CMD>"
  << std::endl << std::endl;
}



void Bustools_sort_Usage() {
  std::cout << "Usage: bustools sort [options] bus-files" << std::endl << std::endl
  << "Options: " << std::endl
  << "-t, --threads         Number of threads to use" << std::endl
  << "-o, --output          File for sorted output" << std::endl
  << "-p, --pipe            Write to standard output" << std::endl
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
  << "-p, --pipe            Write to standard output" << std::endl
  << std::endl;
}

void Bustools_correct_Usage() {
  std::cout << "Usage: bustools correct [options] bus-files" << std::endl << std::endl
  << "Options: " << std::endl
  << "-o, --output          File for corrected bus output" << std::endl
  << "-w, --whitelist       File of whitelisted barcodes to correct to" << std::endl
  << "-p, --pipe            Write to standard output" << std::endl
  << std::endl;
}

void Bustools_count_Usage() {
  std::cout << "Usage: bustools count [options] bus-files" << std::endl << std::endl
  << "Options: " << std::endl
  << "-o, --output          File for corrected bus output" << std::endl
  << "-g, --genemap         File for mapping transcripts to genes" << std::endl
  << "-e, --ecmap           File for mapping equivalence classes to transcripts" << std::endl
  << "-t, --txnames         File with names of transcripts" << std::endl
  << "--genecounts          Aggregate counts to genes only" << std::endl
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

std::vector<int32_t> union_vectors(const std::vector<std::vector<int32_t>> &v) {
  std::vector<int32_t> u;
  for (const auto &vv : v) {
    for (auto t : vv) {
      u.push_back(t);
    }
  }
  
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  return std::move(u);
}

std::vector<int32_t> intersect_vectors(const std::vector<std::vector<int32_t>> &v) {
  std::vector<int32_t> u;

  if (!v.empty()) {
    u = v[0]; // copy
    for (size_t i = 1; i < v.size(); i++) {
      auto &vv = v[i];
      int j = 0;
      int k = 0;
      int l = 0;
      int n = u.size();
      int m = vv.size();
      // u and v are sorted, j,k,l = 0
      while (j < n && l < m) {
        // invariant: u[:k] is the intersection of u[:j] and vv[:l], j <= n, l <= m
        //            u[:j] <= u[j:], vv[:l] <= v[l:], u[j:] is sorted, v[l:] is sorted, u[:k] is sorted
        if (u[j] < vv[l]) {
          j++;
        } else if (u[j] > vv[l]) {
          l++;
        } else {
          // match
          if (k < j) {
            std::swap(u[k], u[j]);
          }
          k++;
          j++;
          l++;
        }
      }
      if (k < n) {
        u.resize(k);
      }    
    }
  }

  return std::move(u);
}

int32_t intersect_ecs(const std::vector<int32_t> &ecs, std::vector<int32_t> &u, std::vector<std::vector<int32_t>> &ecmap, std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> &ecmapinv) {
  if (ecs.empty()) {
    return -1;
  }

  if (ecs[0] < 0 || ecs[0] >= ecmap.size()) {
    return -1;
  }
  
  if (ecs.size() == 1) {
    return ecs[0]; // no work
  }

  u.resize(0);
  auto &v = ecmap[ecs[0]]; // copy
  for (size_t i = 0; i< v.size(); i++) {
    u.push_back(v[i]);
  }

  for (size_t i = 1; i < ecs.size(); i++) {
    if (ecs[i] < 0 || ecs[i] >= ecmap.size()) {
      return -1;
    }
    const auto &v = ecmap[ecs[i]];

    int j = 0;
    int k = 0;
    int l = 0;
    int n = u.size();
    int m = v.size();
    // u and v are sorted, j,k,l = 0
    while (j < n && l < m) {
      // invariant: u[:k] is the intersection of u[:j] and v[:l], j <= n, l <= m
      //            u[:j] <= u[j:], v[:l] <= v[l:], u[j:] is sorted, v[l:] is sorted, u[:k] is sorted
      if (u[j] < v[l]) {
        j++;
      } else if (u[j] > v[l]) {
        l++;
      } else {
        // match
        if (k < j) {
          std::swap(u[k], u[j]);
        }
        k++;
        j++;
        l++;
      }
    }
    if (k < n) {
      u.resize(k);
    }
  }

  if (u.empty()) {
    return -1;
  }
  auto iit = ecmapinv.find(u);
  if (iit == ecmapinv.end()) { 
    // create new equivalence class
    int32_t ec = ecmap.size();
    ecmap.push_back(u);
    ecmapinv.insert({u,ec});
    return ec;
  } else {
    return iit->second;
  }
}


void vt2gene(const std::vector<int32_t> &v, const std::vector<int32_t> &genemap, std::vector<int32_t> &glist) {
  int lastg = -2;
  int n = v.size();

  for (int i = 0; i < n; i++) {
    auto t = v[i];
    auto g = genemap[t];

    if (g != lastg && g != -1) {
      glist.push_back(g);
      lastg = g;
    }
  }

  if (glist.size() > 1) {
    // sort and remove duplicates
    std::sort(glist.begin(), glist.end());
    glist.erase(std::unique(glist.begin(), glist.end()), glist.end());
  }
}



int32_t intersect_ecs_with_genes(const std::vector<int32_t> &ecs, const std::vector<int32_t> &genemap, std::vector<std::vector<int32_t>> &ecmap, std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> &ecmapinv, std::vector<std::vector<int32_t>> &ec2genes, bool assumeIntersectionIsEmpty = true) {
  
  std::vector<std::vector<int32_t>> gu; // per gene transcript results
  std::vector<int32_t> u; // final list of transcripts
  std::vector<int32_t> glist;

  int32_t lastg = -2;
  // todo, replace by intersection of the genelist
  for (const auto ec : ecs) {
    auto g = ec2genes[ec];
    if (g.size() == 1 && g[0] != lastg) {
      glist.push_back(g[0]);
      lastg = g[0];
    } else if (g.size() > 1) {
      lastg = -2;
      for (auto &x : g) {
        glist.push_back(x);
      }
    }
  }
  
  if (glist.empty()) {
    return -1;
  }

  // sort and remove unique
  std::sort(glist.begin(), glist.end());
  glist.erase(std::unique(glist.begin(), glist.end()), glist.end());

  if (glist.size() == 1 && assumeIntersectionIsEmpty) {
    // frequent case, single gene replace with union
    for (auto ec : ecs) {
      for (const auto &t : ecmap[ec]) {      
        u.push_back(t);
      }
    }
    std::sort(u.begin(), u.end());
    u.erase(std::unique(u.begin(), u.end()), u.end());

    // look up ecs based on u
    int32_t ec = -1;
    
    auto it = ecmapinv.find(u);
    if (it != ecmapinv.end()) {
      ec = it->second;              
    } else {
      ec = ecmapinv.size();
      ecmapinv.insert({u,ec});
    }

    return ec; // done
  } else {
    // separate per gene
    for (auto g : glist) {
      gu.clear();
      
      for (auto ec : ecs) {
        std::vector<int32_t> tg;
        for (const auto &t : ecmap[ec]) {
          if (genemap[t] == g) {
            tg.push_back(t);
          }
        }
        if (!tg.empty()) {
          gu.push_back(std::move(tg));
        }
      }
      auto uu = intersect_vectors(gu);

      // if gene intersection is empty, use union
      if (uu.empty()) {
        uu = union_vectors(gu);
      }

      for (auto t : uu) { 
        u.push_back(t);
      }
    }
    int32_t ec = -1;
    auto it = ecmapinv.find(u);
    if (it != ecmapinv.end()) {
      ec = it->second;              
    } else {
      ec = ecmapinv.size();
      ecmapinv.insert({u,ec});
    }
    return ec;
  } 
  
}


void create_ec2genes(const std::vector<std::vector<int32_t>> &ecmap, const std::vector<int32_t> &genemap, std::vector<std::vector<int32_t>> &ec2gene) {
  std::vector<int32_t> u;
  for (int ec = 0; ec < ecmap.size(); ec++) {
    const auto &v = ecmap[ec];    
    vt2gene(v, genemap, u);
    ec2gene.push_back(std::move(u));
  }
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
          while (true) {
            in.read((char*)p, N*sizeof(BUSData));
            size_t rc = in.gcount() / sizeof(BUSData);
            if (rc == 0) {
              break;
            }
            // todo, reserve max memory
            b.insert(b.end(), p, p+rc);
          }
        }
        delete[] p; p = nullptr;
        std::cerr << "Read in " << b.size() << " number of busrecords" << std::endl;

        // todo: replace with radix sort 
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

        std::streambuf *buf = nullptr;
        std::ofstream of;

        if (!opt.stream_out) {
          of.open(opt.output); 
          buf = of.rdbuf();
        } else {
          buf = std::cout.rdbuf();
        }
        std::ostream o(buf);


        char magic[4];      
        uint32_t version = 0;
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
          uint32_t bclen = h.bclen;
          uint32_t umilen = h.umilen;
          int rc = 0;
          while (true) {
            in.read((char*)p, N*sizeof(BUSData));
            size_t rc = in.gcount() / sizeof(BUSData);
            if (rc == 0) {
              break;
            }
            nr += rc;
            for (size_t i = 0; i < rc; i++) {
              o << binaryToString(p[i].barcode, bclen) << "\t" << binaryToString(p[i].UMI,umilen) << "\t" << p[i].ec << "\t" << p[i].count << "\n";        
            }
          }
        }
        delete[] p; p = nullptr;
        if (!opt.stream_out) {
          of.close();
        }
        std::cerr << "Read in " << nr << " number of busrecords" << std::endl;
      } else {
        Bustools_dump_Usage();
        exit(1);
      }
    } else if (cmd == "correct") {      
      if (disp_help) {
        Bustools_correct_Usage();
        exit(0);        
      }
      parse_ProgramOptions_correct(argc-1, argv+1, opt);
      if (check_ProgramOptions_correct(opt)) { //Program options are valid
        
        uint32_t bclen = 0; 
        uint32_t wc_bclen = 0;
        uint32_t umilen = 0;
        BUSHeader h;
        size_t nr = 0;
        size_t N = 100000;
        BUSData* p = new BUSData[N];
        char magic[4];      
        uint32_t version = 0;

        std::ifstream wf(opt.whitelist, std::ios::in);
        std::string line;
        line.reserve(100);
        std::unordered_set<uint64_t> wbc;
        wbc.reserve(100000);
        uint32_t f = 0;
        while(std::getline(wf, line)) {
          if (wc_bclen == 0) {
            wc_bclen = line.size();
          }
          uint64_t bc = stringToBinary(line, f);
          wbc.insert(bc);          
        }
        wf.close();

        std::cerr << "Found " << wbc.size() << " barcodes in the whitelist" << std::endl;

        std::unordered_map<uint64_t, uint64_t> correct;
        // whitelisted barcodes correct to themselves
        for (uint64_t b : wbc) {
          correct.insert({b,b});
        }
        // include hamming distance 1 to all codewords
        std::vector<uint64_t> bad_y;
        for (auto x : wbc) {
          // insert all hamming distance one
          size_t sh = bclen-1;          

          for (size_t i = 0; i < bclen; ++i) {
            for (uint64_t d = 1; d <= 3; d++) {
              uint64_t y = x ^ (d << (2*sh));
              if (correct.find(y) != correct.end()) {
                bad_y.push_back(y);
              } else {
                correct.insert({y,x});
              }
            }                
            sh--;
          }
        }

        // paranoia about error correcting
        int removed_ys = 0;
        for (auto y : bad_y) {
          if (wbc.find(y) == wbc.end()) {
            if (correct.erase(y)>0) { // duplicates are fine
              removed_ys++;
            }
          }
        }

        std::cerr << "Number of hamming dist 1 barcodes = " << correct.size() << std::endl;
        
        std::streambuf *buf = nullptr;
        std::ofstream busf_out;
        
        if (!opt.stream_out) {
          busf_out.open(opt.output , std::ios::out | std::ios::binary);
          buf = busf_out.rdbuf();
        } else {
          buf = std::cout.rdbuf();
        }
        std::ostream bus_out(buf);

        bool outheader_written = false;
        
        nr = 0;
        size_t cc = 0;
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

          if (!outheader_written) {
            writeHeader(bus_out, h);
            outheader_written = true;
          }

          if (bclen == 0) {
            bclen = h.bclen;

            if (bclen != wc_bclen) { 
              std::cerr << "Error: barcode length and whitelist length differ, barcodes = " << bclen << ", whitelist = " << wc_bclen << std::endl
                        << "       check that your whitelist matches the technology used" << std::endl;

              exit(1);
            }
          }
          if (umilen == 0) {
            umilen = h.umilen;
          }

          int rc = 0;
          while (true) {
            in.read((char*)p, N*sizeof(BUSData));
            size_t rc = in.gcount() / sizeof(BUSData);
            if (rc == 0) {
              break;
            }
            nr +=rc;

            for (size_t i = 0; i < rc; i++) {
              bd = p[i];
              auto it = correct.find(bd.barcode);
              if (it != correct.end()) {
                if (bd.barcode != it->second) {
                  bd.barcode = it->second;
                }
                bd.count = 1;
                cc++;
                bus_out.write((char*) &bd, sizeof(bd));
              }
            }
          }
        }

        std::cerr << "Processed " << nr << " bus records, rescued " << cc << " records" << std::endl;
        if (!opt.stream_out) {
          busf_out.close();
        }

        delete[] p; p = nullptr;
      } else {
        Bustools_dump_Usage();
        exit(1);
      }
    } else if (cmd == "count") {
      if (disp_help) {
        Bustools_count_Usage();
        exit(0);        
      }
      parse_ProgramOptions_count(argc-1, argv+1, opt);
      if (check_ProgramOptions_count(opt)) { //Program options are valid
        BUSHeader h;
        size_t nr = 0;
        size_t N = 100000;
        uint32_t bclen = 0;
        BUSData* p = new BUSData[N];

        // read and parse the equivelence class files

        std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> ecmapinv;
        std::vector<std::vector<int32_t>> ecmap;

        std::unordered_map<std::string, int32_t> txnames;
        parseTranscripts(opt.count_txp, txnames);
        std::vector<int32_t> genemap(txnames.size(), -1);
        std::unordered_map<std::string, int32_t> genenames;
        parseGenes(opt.count_genes, txnames, genemap, genenames);
        parseECs(opt.count_ecs, h);
        ecmap = std::move(h.ecs);
        ecmapinv.reserve(ecmap.size());
        for (int32_t ec = 0; ec < ecmap.size(); ec++) {
          ecmapinv.insert({ecmap[ec], ec});
        }
        std::vector<std::vector<int32_t>> ec2genes;        
        create_ec2genes(ecmap, genemap, ec2genes);


        std::ofstream of;
        std::string mtx_ofn = opt.output + ".mtx";
        std::string barcodes_ofn = opt.output + ".barcodes.txt";
        std::string ec_ofn = opt.output + ".ec.txt";
        of.open(mtx_ofn); 

        // write out the initial header
        of << "%%MatrixMarket matrix coordinate real general\n%\n";
        // number of genes
        auto mat_header_pos = of.tellp();
        std::string dummy_header(66, '\n');
        for (int i = 0; i < 33; i++) {
          dummy_header[2*i] = '%';
        }
        of.write(dummy_header.c_str(), dummy_header.size());
    

        size_t n_cols = 0;
        size_t n_rows = 0;
        size_t n_entries = 0;
        std::vector<BUSData> v;
        v.reserve(N);
        uint64_t current_bc = 0xFFFFFFFFFFFFFFFFULL;
        //temporary data
        std::vector<int32_t> ecs;
        ecs.reserve(100);
        std::vector<int32_t> u;
        u.reserve(100);
        std::vector<int32_t> ecs_v;
        ecs_v.reserve(N);
        //barcodes 
        std::vector<uint64_t> barcodes;
        int bad_count = 0;
        int compacted = 0;
        int rescued = 0;

        auto write_barcode_matrix = [&](const std::vector<BUSData> &v) {
          if(v.empty()) {
            return;
          }
          ecs_v.resize(0);
          n_rows+= 1;
          
          barcodes.push_back(v[0].barcode);
          double val = 0.0;
          size_t n = v.size();

          for (size_t i = 0; i < n; ) {
            size_t j = i+1;
            for (; j < n; j++) {
              if (v[i].UMI != v[j].UMI) {
                break;
              }
            }

            // v[i..j-1] share the same UMI
            ecs.resize(0);
            for (size_t k = i; k < j; k++) {
              ecs.push_back(v[k].ec);
            }
            int32_t ec = intersect_ecs(ecs, u, ecmap, ecmapinv);
            if (ec == -1) {
              ec = intersect_ecs_with_genes(ecs, genemap, ecmap, ecmapinv, ec2genes);              
              if (ec == -1) {
                bad_count += j-i;
              } else {
                rescued += j-i;
                ecs_v.push_back(ec);
              }

            } else {
              compacted += j-i-1;
              ecs_v.push_back(ec);
            }
            i = j; // increment
          }
          std::sort(ecs_v.begin(), ecs_v.end());
          size_t m = ecs_v.size();
          for (size_t i = 0; i < m; ) {
            size_t j = i+1;
            for (; j < m; j++) {
              if (ecs_v[i] != ecs_v[j]) {
                break;
              }
            }
            double val = j-i;
            of << n_rows << " " << (ecs_v[i]+1) << " " << val << "\n";
            n_entries++;
            
            i = j; // increment
          }
        };

        for (const auto& infn : opt.files) { 
          std::ifstream inf(infn.c_str(), std::ios::binary);
          parseHeader(inf, h);
          bclen = h.bclen;
          
          int rc = 0;
          while (true) {
            inf.read((char*)p, N*sizeof(BUSData));
            size_t rc = inf.gcount() / sizeof(BUSData);
            nr += rc;
            if (rc == 0) {
              break;
            }

            
            for (size_t i = 0; i < rc; i++) {
              if (p[i].barcode != current_bc) {                 
                // output whatever is in v
                if (!v.empty()) {
                  write_barcode_matrix(v);
                }
                v.clear();
                current_bc = p[i].barcode;
              }
              v.push_back(p[i]);

            }            
          }
          if (!v.empty()) {
            write_barcode_matrix(v);
          }
        }
        delete[] p; p = nullptr;
        n_cols = ecmap.size();
      
        of.close();
        
        std::stringstream ss;
        ss << n_rows << " " << n_cols << " " << n_entries << "\n";
        std::string header = ss.str();
        int hlen = header.size();
        assert(hlen < 66);
        of.open(mtx_ofn, std::ios::binary | std::ios::in | std::ios::out);
        of.seekp(mat_header_pos);
        of.write("%",1);
        of.write(std::string(66-hlen-2,' ').c_str(),66-hlen-2);
        of.write("\n",1);
        of.write(header.c_str(), hlen);
        of.close();

        // write updated ec file
        h.ecs = std::move(ecmap);
        writeECs(ec_ofn, h);
        // write barcode file
        std::ofstream bcof;
        bcof.open(barcodes_ofn);
        for (const auto &x : barcodes) {
          bcof << binaryToString(x, bclen) << "\n";
        }
        bcof.close();
        std::cerr << "bad counts = " << bad_count <<", rescued  =" << rescued << ", compacted = " << compacted << std::endl;

        //std::cerr << "Read in " << nr << " number of busrecords" << std::endl;
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
