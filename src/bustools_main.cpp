#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>
#include <unistd.h>
#include <libgen.h>
#include <string.h>

#ifdef _WIN64
#include <windows.h>
#endif

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

#include "bustools_sort.h"
#include "bustools_count.h"
#include "bustools_whitelist.h"
#include "bustools_project.h"
#include "bustools_inspect.h"
#include "bustools_linker.h"
#include "bustools_capture.h"
#include "bustools_correct.h"
#include "bustools_merge.h"
#include "bustools_extract.h"
#include "bustools_mash.h"
#include "bustools_text.h"
#include "bustools_umicorrect.h"
#include "bustools_predict.h"
#include "bustools_collapse.h"
#include "bustools_clusterhist.h"
#include "bustools_compress.h"
#include "bustools_decompress.h"

int my_mkdir(const char *path, mode_t mode)
{
#ifdef _WIN64
  return mkdir(path);
#else
  return mkdir(path, mode);
#endif
}

bool checkFileExists(const std::string &fn)
{
  struct stat stFileInfo;
  auto intStat = stat(fn.c_str(), &stFileInfo);
  return intStat == 0;
}

bool checkDirectoryExists(const std::string &fn)
{
  struct stat stFileInfo;
  auto intStat = stat(fn.c_str(), &stFileInfo);
  return intStat == 0 && S_ISDIR(stFileInfo.st_mode);
}

bool checkOutputFileValid(const std::string &fn)
{
  std::ofstream of(fn);
  if (of.is_open())
  {
    of.close();
    return checkFileExists(fn);
  }
  else
  {
    return false;
  }
}

std::vector<std::string> parseList(const std::string &s, const std::string &sep = ",")
{
  std::vector<std::string> ret;
  size_t start = 0, end;
  while ((end = s.find(sep, start)) != std::string::npos)
  {
    ret.push_back(s.substr(start, end - start));
    start = end + sep.length();
  }
  end = s.size();
  if (end - start > 0)
  {
    ret.push_back(s.substr(start, end - start));
  }
  return ret;
}

void parse_ProgramOptions_sort(int argc, char **argv, Bustools_opt &opt)
{
  
  const char *opt_string = "t:o:m:T:cuspn";
  
  static struct option long_options[] = {
    {"threads", required_argument, 0, 't'},
    {"output", required_argument, 0, 'o'},
    {"memory", required_argument, 0, 'm'},
    {"temp", required_argument, 0, 'T'},
    {"umi", no_argument, 0, 'u'},
    {"no-flags", no_argument, 0, 'n'},
    {"count", no_argument, 0, 'c'},
    {"flags", no_argument, 0, 'F'},
    {"flags-bc", no_argument, 0, 'f'},
    {"pipe", no_argument, 0, 'p'},
    {0, 0, 0, 0}};
  
  int option_index = 0, c;
  
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1)
  {
    std::string s;
    size_t sh = 0;
    int n = 0;
    switch (c)
    {
      
    case 't':
      opt.threads = atoi(optarg);
      break;
    case 'o':
      opt.output = optarg;
      break;
    case 'n':
      opt.sort_noflag = true; // Not implemented (ideally, would want to have this remove the flag column while sorting)
      break;
    case 'm':
      s = optarg;
      sh = 0;
      n = s.size();
      if (n == 0)
      {
        break;
      }
      switch (s[n - 1])
      {
      case 'm':
      case 'M':
        sh = 20;
        n--;
        break;
      case 'g':
      case 'G':
        sh = 30;
        n--;
        break;
      default:
        sh = 0;
      break;
      }
      opt.max_memory = atoi(s.substr(0, n).c_str());
      opt.max_memory <<= sh;
      break;
    case 'T':
      opt.temp_files = optarg;
      break;
    case 'c':
      opt.type = SORT_COUNT;
      break;
    case 'u':
      opt.type = SORT_UMI;
      break;
    case 'F':
      opt.type = SORT_F;
      break;
    case 'f':
      opt.type = SORT_F_BC;
      break;
    case 'p':
      opt.stream_out = true;
      break;
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  
  // all other arguments are fast[a/q] files to be read
  while (optind < argc)
    opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-")
  {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_merge(int argc, char **argv, Bustools_opt &opt)
{
  const char *opt_string = "o:e:t:";
  static struct option long_options[] = {
    {"output", required_argument, 0, 'o'},
    
    {"ecmap", required_argument, 0, 'e'},
    {"txnames", required_argument, 0, 't'},
    {0, 0, 0, 0}};
  
  int option_index = 0, c;
  
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1)
  {
    
    switch (c)
    {
    case 'o':
      opt.output = optarg;
      break;
    case 't':
      opt.count_txp = optarg;
      break;
    case 'e':
      opt.count_ecs = optarg;
      break;
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  
  while (optind < argc)
    opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-")
  {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_mash(int argc, char **argv, Bustools_opt &opt)
{
  const char *opt_string = "o:";
  
  static struct option long_options[] = {
    {"output", required_argument, 0, 'o'},
    {0, 0, 0, 0}};
  
  int option_index = 0, c;
  
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1)
  {
    
    switch (c)
    {
    case 'o':
      opt.output = optarg;
      break;
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  
  while (optind < argc)
    opt.files.push_back(argv[optind++]);
}

void parse_ProgramOptions_capture(int argc, char **argv, Bustools_opt &opt)
{
  const char *opt_string = "o:xc:e:t:Fsubfp";
  
  static struct option long_options[] = {
    {"output", required_argument, 0, 'o'},
    {"complement", no_argument, 0, 'x'},
    {"capture", required_argument, 0, 'c'},
    {"ecmap", required_argument, 0, 'e'},
    {"txnames", required_argument, 0, 't'},
    {"flags", no_argument, 0, 'F'},
    {"transcripts", no_argument, 0, 's'},
    {"umis", no_argument, 0, 'u'},
    {"barcode", no_argument, 0, 'b'},
    {"combo", no_argument, 0, 'f'},
    {"pipe", no_argument, 0, 'p'},
    {0, 0, 0, 0}};
  
  int option_index = 0, c;
  
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1)
  {
    
    switch (c)
    {
    case 'o':
      opt.output = optarg;
      break;
    case 'x':
      opt.complement = true;
      break;
    case 'c':
      opt.capture = optarg;
      break;
    case 'e':
      opt.count_ecs = optarg;
      break;
    case 't':
      opt.count_txp = optarg;
      break;
    case 'F':
      opt.type = CAPTURE_F;
      break;
    case 's':
      opt.type = CAPTURE_TX;
      break;
    case 'u':
      opt.type = CAPTURE_UMI;
      break;
    case 'b':
      opt.type = CAPTURE_BC;
      break;
    case 'f':
      opt.filter = true;
      break;
    case 'p':
      opt.stream_out = true;
      break;
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  
  while (optind < argc)
    opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-")
  {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_count(int argc, char **argv, Bustools_opt &opt)
{
  const char *opt_string = "o:g:e:t:md:s:";
  int gene_flag = 0;
  int umigene_flag = 0;
  int em_flag = 0;
  int cm_flag = 0;
  int hist_flag = 0;
  int rawcounts_flag = 0;
  static struct option long_options[] = {
    {"output", required_argument, 0, 'o'},
    {"genemap", required_argument, 0, 'g'},
    {"ecmap", required_argument, 0, 'e'},
    {"txnames", required_argument, 0, 't'},
    {"genecounts", no_argument, &gene_flag, 1},
    {"umi-gene", no_argument, &umigene_flag, 1},
    {"multimapping", no_argument, 0, 'm'},
    {"em", no_argument, &em_flag, 1},
    {"cm", no_argument, &cm_flag, 1},
    {"hist", no_argument, &hist_flag, 1},
    {"downsample", required_argument, 0, 'd'},
    {"rawcounts", no_argument, &rawcounts_flag, 1},
    {"split", required_argument, 0, 's'},
    {0, 0, 0, 0}};
  
  int option_index = 0, c;
  
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1)
  {
    
    switch (c)
    {
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
    case 'd':
      opt.count_downsampling_factor = atof(optarg);
      break;
    case 'm':
      opt.count_gene_multimapping = true;
      break;
    case 's':
      opt.count_split = optarg;
      break;
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  if (gene_flag) {
    opt.count_collapse = true;
  }
  if (umigene_flag) {
    opt.umi_gene_collapse = true;
  }
  if (em_flag) {
    opt.count_em = true;
  }
  if (cm_flag) {
    opt.count_cm = true;
  }
  if (hist_flag) {
    opt.count_gen_hist = true;
  }
  if (rawcounts_flag) {
    opt.count_raw_counts = true;
  }
  
  while (optind < argc)
    opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-")
  {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_umicorrect(int argc, char **argv, Bustools_opt& opt) {
  const char* opt_string = "o:g:e:t:p";
  int gene_flag = 0;
  int em_flag = 0;
  int hist_flag = 0;
  static struct option long_options[] = {
    {"output",          required_argument,  0, 'o'},
    {"pipe",            no_argument,        0, 'p'},
    {"genemap",          required_argument,  0, 'g'},
    {"ecmap",          required_argument,  0, 'e'},
    {"txnames",          required_argument,  0, 't'},
    {0,                 0,                  0,  0 }
  };
  
  int option_index = 0, c;
  
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {
    
    //reuse the same flags as for count
    switch (c) {
    case 'o':
      opt.output = optarg;
      break;
    case 'p':
      opt.stream_out = true;
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
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  
  while (optind < argc) opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-") {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_predict(int argc, char **argv, Bustools_opt& opt) {
  const char* opt_string = "o:t:";
  int gene_flag = 0;
  int em_flag = 0;
  int hist_flag = 0;
  static struct option long_options[] = {
    {"output",          required_argument,  0, 'o'},
    {"predict_t",          required_argument,  0, 't'},
    {0,                 0,                  0,  0 }
  };
  
  int option_index = 0, c;
  
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {
    
    switch (c) {
    case 'o':
      opt.output = optarg;
      break;
    case 't':
      opt.predict_t = atof(optarg);
      break;
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  
  if (optind < argc) {
    opt.predict_input = argv[optind++];
  }
}

void parse_ProgramOptions_dump(int argc, char **argv, Bustools_opt &opt)
{
  
  const char *opt_string = "o:pfda";
  
  static struct option long_options[] = {
    {"output", required_argument, 0, 'o'},
    {"pipe", no_argument, 0, 'p'},
    {"flags", no_argument, 0, 'f'},
    {"pad", no_argument, 0, 'd'},
    {"showAll", no_argument, 0, 'a'},
    {0, 0, 0, 0}};
  
  int option_index = 0, c;
  
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1)
  {
    
    switch (c)
    {
    case 'o':
      opt.output = optarg;
      break;
    case 'p':
      opt.stream_out = true;
      break;
    case 'f':
      opt.text_dumpflags = true;
      break;
    case 'd':
      opt.text_dumppad = true;
      break;
    case 'a':
      opt.text_showall = true;
      break;
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  
  // all other arguments are fast[a/q] files to be read
  while (optind < argc)
    opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-")
  {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_fromtext(int argc, char **argv, Bustools_opt& opt) {
  
  const char* opt_string = "o:p";
  
  static struct option long_options[] = {
    {"output",          required_argument,  0, 'o'},
    {"pipe",            no_argument, 0, 'p'},
    {"flags",           no_argument, 0, 'f'},
    {0,                 0,                  0,  0 }
  };
  
  int option_index = 0, c;
  
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) 
  {
    
    switch (c) 
    {
    case 'o':
      opt.output = optarg;
      break;
    case 'p':
      opt.stream_out = true;
      break;
    case '?':
      opt.parse_error = true;
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

void parse_ProgramOptions_correct(int argc, char **argv, Bustools_opt &opt)
{
 
  int nocorrect_flag = 0; 
  const char *opt_string = "o:w:d:spr";
  static struct option long_options[] = {
    {"output", required_argument, 0, 'o'},
    {"whitelist", required_argument, 0, 'w'},
    {"onlist", required_argument, 0, 'w'},
    {"dump", required_argument, 0, 'd'},
    {"split", no_argument, 0, 's'},
    {"pipe", no_argument, 0, 'p'},
    {"replace", no_argument, 0, 'r'},
    {"nocorrect", no_argument, &nocorrect_flag, 1},
    {0, 0, 0, 0}};
  
  int option_index = 0, c;
  
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1)
  {
    
    switch (c)
    {
    case 'o':
      opt.output = optarg;
      break;
    case 'w':
      opt.whitelist = optarg;
      break;
    case 'd':
      opt.dump = optarg;
      opt.dump_bool = true;
      break;
    case 's':
      opt.split_correct = true;
      break;
    case 'p':
      opt.stream_out = true;
      break;
    case 'r':
      opt.barcode_replacement = true;
      break;
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  
  //hard code options for now
  opt.ec_d = 1;
  opt.ec_dmin = 3;
  
  // all other arguments are fast[a/q] files to be read
  while (optind < argc)
    opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-")
  {
    opt.stream_in = true;
  }
  if (nocorrect_flag) {
    opt.no_correct = true;
  }
}

void parse_ProgramOptions_whitelist(int argc, char **argv, Bustools_opt &opt)
{
  
  /* Parse options. */
  const char *opt_string = "o:f:h";
  
  static struct option long_options[] = {
    {"output", required_argument, 0, 'o'},
    {"threshold", required_argument, 0, 'f'},
    {0, 0, 0, 0}};
  
  int option_index = 0, c;
  
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1)
  {
    switch (c)
    {
    case 'o':
      opt.output = optarg;
      break;
    case 'f':
      opt.threshold = atoi(optarg);
      break;
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  
  /* All other arguments are (sorted) BUS files. */
  while (optind < argc)
    opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-")
  {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_project(int argc, char **argv, Bustools_opt &opt)
{
  
  /* Parse options. */
  const char *opt_string = "o:m:e:t:s:Fbup";
  
  static struct option long_options[] = {
    {"output", required_argument, 0, 'o'},
    {"map", required_argument, 0, 'm'},
    {"ecmap", required_argument, 0, 'e'},
    {"txnames", required_argument, 0, 't'},
    {"flags", no_argument, 0, 'F'},
    {"barcode", no_argument, 0, 'b'},
    {"umi", no_argument, 0, 'u'},
    {"transcripts", optional_argument, 0, 's'},
    {"pipe", no_argument, 0, 'p'},
    {0, 0, 0, 0}};
  
  int option_index = 0, c;
  
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1)
  {
    switch (c)
    {
    case 'o':
      opt.output = optarg;
      break;
    case 'm':
      opt.map = optarg;
      break;
    case 'e':
      opt.count_ecs = optarg;
      break;
    case 't':
      opt.count_txp = optarg;
      break;
    case 'F':
      opt.type = PROJECT_F;
      break;
    case 'b':
      opt.type = PROJECT_BC;
      break;
    case 'u':
      opt.type = PROJECT_UMI;
      break;
    case 's':
      opt.type = PROJECT_TX;
      opt.output_folder = optarg;
      break;
    case 'p':
      opt.stream_out = true;
      break;
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  
  /* All other arguments are (sorted) BUS files. */
  while (optind < argc)
    opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-")
  {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_inspect(int argc, char **argv, Bustools_opt &opt)
{
  
  /* Parse options. */
  const char *opt_string = "o:e:w:p";
  
  static struct option long_options[] = {
    {"output", required_argument, 0, 'o'},
    {"ecmap", required_argument, 0, 'e'},
    {"whitelist", required_argument, 0, 'w'},
    {"onlist", required_argument, 0, 'w'},
    {"pipe", no_argument, 0, 'p'},
    {0, 0, 0, 0}};
  
  int option_index = 0, c;
  
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1)
  {
    switch (c)
    {
    case 'o':
      opt.output = optarg;
      break;
    case 'e':
      opt.count_ecs = optarg;
      break;
    case 'w':
      opt.whitelist = optarg;
      break;
    case 'p':
      opt.stream_out = true;
      break;
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  
  /* All other arguments are (sorted) BUS files. */
  while (optind < argc)
    opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-")
  {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_linker(int argc, char **argv, Bustools_opt &opt)
{
  
  /* Parse options. */
  const char *opt_string = "o:s:e:p";
  
  static struct option long_options[] = {
    {"output", required_argument, 0, 'o'},
    {"start", required_argument, 0, 's'},
    {"end", required_argument, 0, 'e'},
    {"pipe", no_argument, 0, 'p'},
    {0, 0, 0, 0}};
  
  int option_index = 0, c;
  
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1)
  {
    switch (c)
    {
    case 'o':
      opt.output = optarg;
      break;
    case 's':
      opt.start = std::stoi(optarg);
      break;
    case 'e':
      opt.end = std::stoi(optarg);
      break;
    case 'p':
      opt.stream_out = true;
      break;
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  
  /* All other arguments are (sorted) BUS files. */
  while (optind < argc)
    opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-")
  {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_collapse(int argc, char** argv, Bustools_opt& opt) {
  const char* opt_string = "o:pg:e:t:";
  int gene_flag = 0;
  int em_flag = 0;
  static struct option long_options[] = {
    {"output",          required_argument,  0, 'o'},
    {"pipe",            no_argument, 0, 'p'},
    {"genemap",         required_argument,  0, 'g'},
    {"ecmap",           required_argument,  0, 'e'},
    {"txnames",         required_argument,  0, 't'},
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
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  
  while (optind < argc) opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-") {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_clusterhist(int argc, char** argv, Bustools_opt& opt) {
  const char* opt_string = "o:pg:e:t:c:";
  int gene_flag = 0;
  int em_flag = 0;
  static struct option long_options[] = {
    {"output",          required_argument,  0, 'o'},
    {"pipe",            no_argument, 0, 'p'},
    {"genemap",         required_argument,  0, 'g'},
    {"ecmap",           required_argument,  0, 'e'},
    {"txnames",         required_argument,  0, 't'},
    {"clusterfile",     required_argument,  0, 'c'},
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
    case 'c':
      opt.cluster_input_file = optarg;
      break;
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  
  while (optind < argc) opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-") {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_extract(int argc, char **argv, Bustools_opt &opt)
{
  
  /* Parse options. */
  const char *opt_string = "o:f:N:p";
  
  static struct option long_options[] = {
    {"output", required_argument, 0, 'o'},
    {"fastq", required_argument, 0, 'f'},
    {"nFastqs", required_argument, 0, 'N'},
    {"pipe", no_argument, 0, 'p'},
    {"exclude", no_argument, 0, 'x'},
    {"include", no_argument, 0, 'i'},
    {0, 0, 0, 0}};
  
  int option_index = 0, c;
  
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1)
  {
    switch (c)
    {
    case 'o':
      opt.output = optarg;
      break;
    case 'f':
      opt.fastq = parseList(optarg);
      break;
    case 'N':
      opt.nFastqs = std::stoi(optarg);
      break;
    case 'p':
      opt.stream_out = true;
      break;
    case '?':
      opt.parse_error = true;
      break;
    case 'x':
      opt.extract_exclude = true;
      opt.extract_include = false;
      break;
    case 'i':
      opt.extract_include = true;
      break;
    default:
      break;
    }
  }
  
  /* All other arguments are (sorted) BUS files. */
  while (optind < argc)
    opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-")
  {
    opt.stream_in = true;
  }
}

/**
 * @brief Parse command line arguments for bustools inflate.
 *
 * @param argc
 * @param argv
 * @param opt
 * @return bool true iff requested from the command line.
 */
bool parse_ProgramOptions_inflate(int argc, char **argv, Bustools_opt &opt)
{
  const char *opt_string = "o:ph";
  bool print_usage = false;
  static struct option long_options[] = {
      {"output", required_argument, 0, 'o'},
      {"pipe", no_argument, 0, 'p'},
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0},
  };
  int option_index = 0, c;
  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1)
  {
    switch (c)
    {
    case 'o':
      opt.output = optarg;
      break;
    case 'p':
      opt.stream_out = true;
      break;
    case 'h':
      print_usage = true;
      break;
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }

  while (optind < argc)
    opt.files.push_back(argv[optind++]);

  if (opt.files.size() == 1 && opt.files[0] == "-")
  {
    opt.stream_in = true;
  }
  return print_usage;
}

/**
 * @brief Parse command line arguments for bustools compress.
 *
 * @param argc
 * @param argv
 * @param opt
 * @return bool true iff requested from the command line.
 */
bool parse_ProgramOptions_compress(int argc, char **argv, Bustools_opt &opt)
{
  const char *opt_string = "N:Lo:pP:h:i::";
  const bool lossy_umi_enabled = false;

  static struct option long_options[] = {
      {"chunk-size", required_argument, 0, 'N'},
      {"lossy-umi", no_argument, 0, 'L'},
      {"output", required_argument, 0, 'o'},
      {"pipe", no_argument, 0, 'p'},
      {"pfd-size", required_argument, 0, 'P'},
      {"index", optional_argument, 0, 'i'},
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}};
  int option_index = 0, c;
  bool print_usage = false;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1)
  {
    switch (c)
    {
    case 'i':
    {
      if (optarg == NULL && optind < argc - 1 && argv[optind][0] != '-')
      {
      optarg = argv[optind++];
      }
      if (optarg == NULL)
      {
        opt.busz_index.assign("output.busz.idx");
      }
      else
      {
        opt.busz_index.assign(optarg);
      }
      break;
    }
    case 'N':
      opt.chunk_size = atoi(optarg);
      break;
    case 'L':
      if (!lossy_umi_enabled)
      std::cerr << "Lossy UMI not yet implemented. Using lossless instead." << std::endl;
      opt.lossy_umi = lossy_umi_enabled;
      break;
    case 'o':
      opt.output = optarg;
      break;
    case 'p':
      opt.stream_out = true;
      break;
    case 'P':
      opt.pfd_blocksize = atoi(optarg);
      break;
    case 'h':
      print_usage = true;
      break;
    case '?':
      opt.parse_error = true;
      break;
    default:
      break;
    }
  }
  // all other arguments are fast[a/q] files to be read
  while (optind < argc)
    opt.files.push_back(argv[optind++]);

  if (opt.files.size() == 1 && opt.files[0] == "-")
  {
    opt.stream_in = true;
  }
  return print_usage;
}

bool check_ProgramOptions_sort(Bustools_opt &opt)
{
  
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  size_t max_threads = std::thread::hardware_concurrency();
  
  if (opt.threads <= 0)
  {
    std::cerr << "Error: Number of threads cannot be less than or equal to 0" << std::endl;
    ret = false;
  }
  else if (opt.threads > max_threads)
  {
    std::cerr << "Warning: Number of threads cannot be greater than or equal to " << max_threads
              << ". Setting number of threads to " << max_threads << std::endl;
    opt.threads = max_threads;
  }
  
  if (!opt.stream_out)
  {
    if (opt.output.empty())
    {
      std::cerr << "Error: missing output file" << std::endl;
      ret = false;
    }
    else if (!checkOutputFileValid(opt.output))
    {
      std::cerr << "Error: unable to open output file" << std::endl;
      ret = false;
    }
  }
  
  if (opt.max_memory < 1ULL << 26)
  {
    if (opt.max_memory < 128)
    {
      std::cerr << "Warning: low number supplied for maximum memory usage without M or G suffix\n  interpreting this as " << opt.max_memory << "Gb" << std::endl;
      opt.max_memory <<= 30;
    }
    else
    {
      std::cerr << "Warning: low number supplied for maximum memory, defaulting to 64Mb" << std::endl;
      opt.max_memory = 1ULL << 26; // 64Mb is absolute minimum
    }
  }
  
  if (opt.temp_files.empty())
  {
    if (opt.stream_out)
    {
      std::cerr << "Error: temporary location -T must be set when using streaming output" << std::endl;
      ret = false;
    }
    else
    {
      opt.temp_files = opt.output + ".";
    }
  }
  else
  {
    
    if (checkDirectoryExists(opt.temp_files))
    {
      // if it is a directory, create random file prefix
      opt.temp_files += "/bus.sort." + std::to_string(getpid()) + ".";
    }
    else
    {
      if (opt.temp_files.at(opt.temp_files.size() - 1) == '/')
      {
        // try creating the directory
        if (my_mkdir(opt.temp_files.c_str(), 0777) == 0)
        {
          //
          opt.temp_files += "/bus.sort." + std::to_string(getpid()) + ".";
        }
        else
        {
          std::cerr << "Error: directory " << opt.temp_files << " does not exist and could not be created. Check that the parent directory exists and you have write permissions." << std::endl;
          ret = false;
        }
      }
      else
      {
        int n = opt.temp_files.size();
        if (opt.temp_files[n - 1] != '.')
        {
          opt.temp_files += '.';
        }
        // assume the directory exist
      }
    }
  }
  
  if (ret)
  {
    FILE *tmpf = fopen(opt.temp_files.c_str(), "a+");
    if (tmpf == nullptr)
    {
      char *dirn = dirname(strdup(opt.temp_files.c_str()));
      std::cerr << "Error: Could not create a temporary file in directory " << dirn << ". check that the directory exists and if you have permissions to write in it" << std::endl;
      free(dirn);
      ret = false;
    }
    remove(opt.temp_files.c_str());
  }
  
  if (opt.files.size() == 0)
  {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  }
  else if (!opt.stream_in)
  {
    for (const auto &it : opt.files)
    {
      if (!checkFileExists(it))
      {
        std::cerr << "Error: File not found, " << it << std::endl;
        ret = false;
      }
    }
  }
  
  return ret;
}

bool check_ProgramOptions_merge(Bustools_opt &opt)
{
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  // check for output directory
  if (!opt.stream_out)
  {
    if (opt.output.empty())
    {
      std::cerr << "Error: missing output directory" << std::endl;
      ret = false;
    }
    else
    {
      // check if output directory exists or if we can create it
      struct stat stFileInfo;
      auto intStat = stat(opt.output.c_str(), &stFileInfo);
      if (intStat == 0)
      {
        // file/dir exits
        if (!S_ISDIR(stFileInfo.st_mode))
        {
          std::cerr << "Error: file " << opt.output << " exists and is not a directory" << std::endl;
          ret = false;
        }
      }
      else
      {
        // create directory
        if (my_mkdir(opt.output.c_str(), 0777) == -1)
        {
          std::cerr << "Error: could not create directory " << opt.output << std::endl;
          ret = false;
        }
      }
    }
  }
  
  if (opt.files.size() == 0)
  {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  }
  else
  {
    if (!opt.stream_in)
    {
      for (const auto &it : opt.files)
      {
        if (!checkFileExists(it))
        {
          std::cerr << "Error: File not found, " << it << std::endl;
          ret = false;
        }
      }
    }
  }
  
  if (opt.count_ecs.size() == 0)
  {
    std::cerr << "Error: missing equivalence class mapping file" << std::endl;
    ret = false;
  }
  else
  {
    if (!checkFileExists(opt.count_ecs))
    {
      std::cerr << "Error: File not found " << opt.count_ecs << std::endl;
      ret = false;
    }
  }
  
  if (opt.count_txp.size() == 0)
  {
    std::cerr << "Error: missing transcript name file" << std::endl;
    ret = false;
  }
  else
  {
    if (!checkFileExists(opt.count_txp))
    {
      std::cerr << "Error: File not found " << opt.count_txp << std::endl;
      ret = false;
    }
  }
  
  return ret;
}

bool check_ProgramOptions_mash(Bustools_opt &opt)
{
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  if (opt.output.empty())
  {
    std::cerr << "Error: missing output directory" << std::endl;
  }
  else
  {
    // check if output directory exists or if we can create it
    struct stat stFileInfo;
    auto intStat = stat(opt.output.c_str(), &stFileInfo);
    if (intStat == 0)
    {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode))
      {
        std::cerr << "Error: file " << opt.output << " exists and is not a directory" << std::endl;
        ret = false;
      }
    }
    else
    {
      // create directory
      if (my_mkdir(opt.output.c_str(), 0777) == -1)
      {
        std::cerr << "Error: could not create directory " << opt.output << std::endl;
        ret = false;
      }
    }
  }
  
  if (opt.files.size() == 0)
  {
    std::cerr << "Error: Missing BUS input directories" << std::endl;
    ret = false;
  }
  else
  {
    for (const auto &it : opt.files)
    {
      if (!checkDirectoryExists(it))
      {
        std::cerr << "Error: directory " << it << " does not exist or is not a directory" << std::endl;
        ret = false;
      }
    }
    if (ret)
    {
      // check for output.bus and matrix.ec in each of the directories
      for (const auto &it : opt.files)
      {
        if (!checkFileExists(it + "/output.bus"))
        {
          std::cerr << "Error: file " << it << "/output.bus not found" << std::endl;
        }
        
        if (!checkFileExists(it + "/matrix.ec"))
        {
          std::cerr << "Error: file " << it << "/matrix.ec not found" << std::endl;
        }
      }
    }
  }
  
  return ret;
}

bool check_ProgramOptions_dump(Bustools_opt &opt)
{
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  if (!opt.stream_out)
  {
    if (opt.output.empty())
    {
      std::cerr << "Error: missing output file" << std::endl;
      ret = false;
    }
    else if (!checkOutputFileValid(opt.output))
    {
      std::cerr << "Error: unable to open output file" << std::endl;
      ret = false;
    }
  }
  
  if (opt.files.size() == 0)
  {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  }
  else if (!opt.stream_in)
  {
    for (const auto &it : opt.files)
    {
      if (!checkFileExists(it))
      {
        std::cerr << "Error: File not found, " << it << std::endl;
        ret = false;
      }
    }
  }
  
  return ret;
}

bool check_ProgramOptions_fromtext(Bustools_opt& opt) 
{
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  if (!opt.stream_out) 
  {
    if (opt.output.empty()) 
    {
      std::cerr << "Error: missing output file" << std::endl;
      ret = false;
    } 
    else if (!checkOutputFileValid(opt.output)) 
    {
      std::cerr << "Error: unable to open output file" << std::endl;
      ret = false;
    }
  } 
  
  
  if (opt.files.size() == 0) 
  {
    std::cerr << "Error: Missing input files" << std::endl;
    ret = false;
  } 
  else if (!opt.stream_in) 
  {    
    for (const auto& it : opt.files) 
    {  
      if (!checkFileExists(it)) 
      {
        std::cerr << "Error: File not found, " << it << std::endl;
        ret = false;
      }
    }
  }
  
  return ret;
}

bool check_ProgramOptions_capture(Bustools_opt &opt)
{
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  if (opt.filter)
  {
    // check if output directory exists or if we can create it
    if (!checkDirectoryExists(opt.output))
    {
      std::cerr << "Error: file " << opt.output << " exists and is not a directory" << std::endl;
      ret = false;
    }
    else
    {
      // create directory
      if (my_mkdir(opt.output.c_str(), 0777) == -1)
      {
        std::cerr << "Error: could not create directory " << opt.output << std::endl;
        ret = false;
      }
    }
  }
  else if (!opt.stream_out)
  {
    if (opt.output.empty())
    {
      std::cerr << "Error: missing output file" << std::endl;
      ret = false;
    }
    else if (!checkOutputFileValid(opt.output))
    {
      std::cerr << "Error: unable to open output file" << std::endl;
      ret = false;
    }
  }
  
  if (opt.capture.empty())
  {
    std::cerr << "Error: missing capture list" << std::endl;
    ret = false;
  }
  else
  {
    if (!checkFileExists(opt.capture))
    {
      std::cerr << "Error: File not found, " << opt.capture << std::endl;
      ret = false;
    }
  }
  
  if (opt.type == CAPTURE_NONE)
  {
    std::cerr << "Error: capture list type must be specified (one of -s, -u, or -b)" << std::endl;
    ret = false;
  }
  
  if (opt.type == CAPTURE_TX)
  {
    if (opt.count_ecs.size() == 0)
    {
      std::cerr << "Error: missing equialence class mapping file" << std::endl;
    }
    else
    {
      if (!checkFileExists(opt.count_ecs))
      {
        std::cerr << "Error: File not found " << opt.count_ecs << std::endl;
        ret = false;
      }
    }
    
    if (opt.count_txp.size() == 0)
    {
      std::cerr << "Error: missing transcript name file" << std::endl;
    }
    else
    {
      if (!checkFileExists(opt.count_txp))
      {
        std::cerr << "Error: File not found " << opt.count_txp << std::endl;
        ret = false;
      }
    }
  }
  
  if (opt.files.size() == 0)
  {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  }
  else if (!opt.stream_in)
  {
    for (const auto &it : opt.files)
    {
      if (!checkFileExists(it))
      {
        std::cerr << "Error: File not found, " << it << std::endl;
        ret = false;
      }
    }
  }
  
  if (opt.filter && (opt.complement || opt.type != CAPTURE_TX))
  {
    std::cerr << "Warning: filter only meaningful without complement flag, and to"
              << " capture transcripts; no new ec file will be generated" << std::endl;
    opt.filter = false;
  }
  
  return ret;
}

bool check_ProgramOptions_correct(Bustools_opt &opt)
{
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  if (!opt.stream_out)
  {
    if (opt.output.empty())
    {
      std::cerr << "Error: missing output file" << std::endl;
      ret = false;
    }
    else if (!checkOutputFileValid(opt.output))
    {
      std::cerr << "Error: unable to open output file" << std::endl;
      ret = false;
    }
  }
  
  if (opt.files.size() == 0)
  {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  }
  else
  {
    if (!opt.stream_in)
    {
      for (const auto &it : opt.files)
      {
        if (!checkFileExists(it))
        {
          std::cerr << "Error: File not found, " << it << std::endl;
          ret = false;
        }
      }
    }
  }
  
  if (opt.whitelist.size() == 0)
  {
    std::cerr << "Error: Missing on-list file" << std::endl;
    ret = false;
  }
  else
  {
    if (!checkFileExists(opt.whitelist))
    {
      std::cerr << "Error: File not found " << opt.whitelist << std::endl;
      ret = false;
    }
  }
  if (opt.dump_bool)
  {
    if (opt.dump.empty())
    {
      std::cerr << "Error: dump file not specified" << std::endl;
      ret = false;
    }
  }
  
  return ret;
}

bool check_ProgramOptions_count(Bustools_opt &opt)
{
  bool ret = true;

  if (opt.parse_error) {
      ret = false;
  }
  
  // check for output directory
  if (opt.output.empty()) {
    std::cerr << "Error: Missing output directory" << std::endl;
    ret = false;
  }
  else {
    bool isDir = false;
    if (checkDirectoryExists(opt.output)) {
      isDir = true;
    }
    else {
      if (opt.output.at(opt.output.size() - 1) == '/') {
        if (my_mkdir(opt.output.c_str(), 0777) == -1) {
          std::cerr << "Error: could not create directory " << opt.output << std::endl;
          ret = false;
        }
        else {
          isDir = true;
        }
      }
    }
    
    if (isDir) {
      opt.output += "output";
    }
  }
  
  if (opt.count_em && opt.count_gene_multimapping) {
    std::cerr << "Error: EM algorithm and counting multimapping reads are incompatible" << std::endl;
    ret = false;
  }
  
  if (opt.count_em && opt.count_cm) {
    std::cerr << "Error: EM algorithm and counting multiplicites are incompatible" << std::endl;
    ret = false;
  }
  
  if (opt.umi_gene_collapse && opt.count_cm) {
    std::cerr << "Error: Gene-level collapsing of UMIs and counting multiplicites are incompatible" << std::endl;
    ret = false;
  }
  
  if (opt.umi_gene_collapse && (opt.count_raw_counts || opt.count_gen_hist || opt.count_downsampling_factor != 1.0)) {
    std::cerr << "Error: Gene-level collapsing of UMIs is currently incompatible with --hist, --downsample, or --rawcounts" << std::endl;
    ret = false;
  }
  
  if (opt.count_cm && (opt.count_raw_counts || opt.count_gen_hist || opt.count_downsampling_factor != 1.0)) {
    std::cerr << "Error: Counting multiplicites is incompatible with --hist, --downsample, or --rawcounts" << std::endl;
    ret = false;
  }
  
  if (opt.count_raw_counts && opt.count_em)  {
    std::cerr << "Error: Counting raw counts are not supported for the EM algorithm" << std::endl;
    ret = false;
  }
  
  if (opt.count_raw_counts && !opt.count_collapse)  {
    std::cerr << "Error: Raw counts are currently only supported for gene counting, not ec counting." << std::endl;
    ret = false;
  }
  
  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  }
  else {
    if (!opt.stream_in) {
      for (const auto &it : opt.files) {
        if (!checkFileExists(it)) {
          std::cerr << "Error: File not found, " << it << std::endl;
          ret = false;
        }
      }
    }
  }
  
  if (opt.count_genes.size() == 0) {
    std::cerr << "Error: missing gene mapping file" << std::endl;
    ret = false;
  }
  else {
    if (!checkFileExists(opt.count_genes))
    {
      std::cerr << "Error: File not found " << opt.count_genes << std::endl;
      ret = false;
    }
  }
  
  if (opt.count_ecs.size() == 0) {
    std::cerr << "Error: missing equivalence class mapping file" << std::endl;
    ret = false;
  }
  else {
    if (!checkFileExists(opt.count_ecs))
    {
      std::cerr << "Error: File not found " << opt.count_ecs << std::endl;
      ret = false;
    }
  }
  
  if (opt.count_txp.size() == 0) {
    std::cerr << "Error: missing transcript name file" << std::endl;
    ret = false;
  }
  else {
    if (!checkFileExists(opt.count_txp)) {
      std::cerr << "Error: File not found " << opt.count_txp << std::endl;
      ret = false;
    }
  }
  
  if (opt.count_split.size() != 0) {
    if (!checkFileExists(opt.count_split)) {
      std::cerr << "Error: File not found " << opt.count_split << std::endl;
      ret = false;
    }
    if (opt.count_em) {
      std::cerr << "Cannot use -s with --em" << std::endl;
      ret = false;
    }
  }
  
  return ret;
}

bool check_ProgramOptions_predict(Bustools_opt& opt) {
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  // check for output directory
  if (opt.output.empty()) {
    std::cerr << "Error: Missing output directory" << std::endl;
    ret = false;
  } else {
    bool isDir = false;
    if (checkDirectoryExists(opt.output)) {
      isDir = true;
    } else {
      if (opt.output.at(opt.output.size()-1) == '/') {
        if (my_mkdir(opt.output.c_str(), 0777) == -1) {
          std::cerr << "Error: could not create directory " << opt.output << std::endl;
          ret = false;
        } else {
          isDir = true;
        }
      }
    }
    
    if (isDir) {
      opt.output += "output";
    }
  }
  
  // check for input directory
  if (opt.output.empty()) {
    std::cerr << "Error: Missing input" << std::endl;
    ret = false;
  } else {
    bool isDir = false;
    if (checkDirectoryExists(opt.predict_input)) {
      isDir = true;
    } else {
      if (opt.predict_input.at(opt.predict_input.size()-1) == '/') {
        isDir = true;
      }
    }
    
    if (isDir) {
      opt.predict_input += "output"; //to match the output from count
    }
    
    if ( !checkFileExists(opt.predict_input + ".mtx") ) {
      std::cerr << "Error: Matrix file missing: " << opt.predict_input + ".mtx" << std::endl;
      ret = false;
    }
    if ( !checkFileExists(opt.predict_input + ".genes.txt") ) {
      std::cerr << "Error: Genes file missing: " << opt.predict_input + ".genes.txt" << std::endl;
      ret = false;
    }
    if ( !checkFileExists(opt.predict_input + ".barcodes.txt") ) {
      std::cerr << "Error: Barcodes file missing: " << opt.predict_input + ".barcodes.txt" << std::endl;
      ret = false;
    }
    if ( !checkFileExists(opt.predict_input + ".hist.txt") ) {
      std::cerr << "Error: CPU histograms file missing: " << opt.predict_input + ".hist.txt. Did you forget the --hist flag when running count?" << std::endl;
      ret = false;
    }
  }
  
  if (opt.predict_t == 0.0) {
    std::cerr << "Error: Prediction range not set." << std::endl;
    ret = false;
  }
  
  return ret;
}

bool check_ProgramOptions_umicorrect(Bustools_opt& opt) {
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  if (!opt.stream_out) {
    if (opt.output.empty()) {
      std::cerr << "Error: missing output file" << std::endl;
      ret = false;
    } else if (!checkOutputFileValid(opt.output)) {
      std::cerr << "Error: unable to open output file" << std::endl;
      ret = false;
    }
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
    ret = false;
  } else {
    if (!checkFileExists(opt.count_genes)) {
      std::cerr << "Error: File not found " << opt.count_genes << std::endl;
      ret = false;
    }
  }
  
  if (opt.count_ecs.size() == 0) {
    std::cerr << "Error: missing equivalence class mapping file" << std::endl;
    ret = false;
  } else {
    if (!checkFileExists(opt.count_ecs)) {
      std::cerr << "Error: File not found " << opt.count_ecs << std::endl;
      ret = false;
    }
  }
  
  if (opt.count_txp.size() == 0) {
    std::cerr << "Error: missing transcript name file" << std::endl;
    ret = false;
  } else {
    if (!checkFileExists(opt.count_txp)) {
      std::cerr << "Error: File not found " << opt.count_txp << std::endl;
      ret = false;
    }
  }
  
  return ret;
}

bool check_ProgramOptions_whitelist(Bustools_opt &opt)
{
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  if (opt.output.empty())
  {
    std::cerr << "Error: missing output file" << std::endl;
    ret = false;
  }
  else if (!checkOutputFileValid(opt.output))
  {
    std::cerr << "Error: unable to open output file" << std::endl;
    ret = false;
  }
  
  if (opt.files.size() == 0)
  {
    std::cerr << "Error: Missing BUS input file" << std::endl;
    ret = false;
  }
  else if (opt.files.size() == 1)
  {
    if (!opt.stream_in)
    {
      for (const auto &it : opt.files)
      {
        if (!checkFileExists(it))
        {
          std::cerr << "Error: File not found, " << it << std::endl;
          ret = false;
        }
      }
    }
  }
  else
  {
    std::cerr << "Error: Only one input file allowed" << std::endl;
    ret = false;
  }
  
  if (opt.threshold < 0)
  { // threshold = 0 for no threshold
    std::cerr << "Error: Threshold cannot be less than or equal to 0" << std::endl;
    ret = false;
  }
  
  return ret;
}

bool check_ProgramOptions_project(Bustools_opt &opt)
{
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  if (opt.output.empty())
  {
    std::cerr << "Error: Missing output file" << std::endl;
  }
  else if (!checkOutputFileValid(opt.output))
  {
    std::cerr << "Error: unable to open output file" << std::endl;
    ret = false;
  }
  if (opt.type == PROJECT_TX)
  {
    if (opt.output_folder.empty())
    {
      std::cerr << "Error: Missing output folder" << std::endl;
    }
    else
    {
      // check if output directory exists or if we can create it
      struct stat stFileInfo;
      auto intStat = stat(opt.output_folder.c_str(), &stFileInfo);
      if (intStat == 0)
      {
        // file/dir exits
        if (!S_ISDIR(stFileInfo.st_mode))
        {
          std::cerr << "Error: file " << opt.output_folder << " exists and is not a directory" << std::endl;
          ret = false;
        }
      }
      else
      {
        // create directory
        if (my_mkdir(opt.output_folder.c_str(), 0777) == -1)
        {
          std::cerr << "Error: could not create directory " << opt.output_folder << std::endl;
          ret = false;
        }
      }
    }
  }
  
  if (opt.files.size() == 0)
  {
    std::cerr << "Error: Missing BUS input file" << std::endl;
    ret = false;
  }
  else if (opt.files.size() == 1)
  {
    if (!opt.stream_in)
    {
      for (const auto &it : opt.files)
      {
        if (!checkFileExists(it))
        {
          std::cerr << "Error: File not found, " << it << std::endl;
          ret = false;
        }
      }
    }
  }
  else
  {
    std::cerr << "Error: Only one input file allowed" << std::endl;
    ret = false;
  }
  
  if (opt.map.size() == 0)
  {
    std::cerr << "Error: missing mapping file" << std::endl;
  }
  else
  {
    if (!checkFileExists(opt.map))
    {
      std::cerr << "Error: File not found " << opt.map << std::endl;
      ret = false;
    }
  }
  
  if (opt.count_ecs.size() == 0)
  {
    std::cerr << "Error: missing equivalence class mapping file" << std::endl;
  }
  else
  {
    if (!checkFileExists(opt.count_ecs))
    {
      std::cerr << "Error: File not found " << opt.count_ecs << std::endl;
      ret = false;
    }
  }
  
  if (opt.count_txp.size() == 0)
  {
    std::cerr << "Error: missing transcript name file" << std::endl;
  }
  else
  {
    if (!checkFileExists(opt.map))
    {
      std::cerr << "Error: File not found " << opt.count_txp << std::endl;
      ret = false;
    }
  }
  
  return ret;
}

bool check_ProgramOptions_inspect(Bustools_opt &opt)
{
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  if (opt.output.size() && !checkOutputFileValid(opt.output))
  {
    std::cerr << "Error: unable to open output file" << std::endl;
    ret = false;
  }
  
  if (opt.files.size() == 0)
  {
    std::cerr << "Error: Missing BUS input file" << std::endl;
    ret = false;
  }
  else if (opt.files.size() == 1)
  {
    if (!opt.stream_in)
    {
      for (const auto &it : opt.files)
      {
        if (!checkFileExists(it))
        {
          std::cerr << "Error: File not found, " << it << std::endl;
          ret = false;
        }
      }
    }
  }
  else
  {
    std::cerr << "Error: Only one input file allowed" << std::endl;
    ret = false;
  }
  
  if (opt.count_ecs.size())
  {
    if (!checkFileExists(opt.count_ecs))
    {
      std::cerr << "Error: File not found " << opt.count_ecs << std::endl;
      ret = false;
    }
  }
  
  if (opt.whitelist.size())
  {
    if (!checkFileExists(opt.whitelist))
    {
      std::cerr << "Error: File not found " << opt.whitelist << std::endl;
      ret = false;
    }
  }
  
  return ret;
}

bool check_ProgramOptions_linker(Bustools_opt &opt)
{
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  if (!opt.stream_out)
  {
    if (opt.output.empty())
    {
      std::cerr << "Error: missing output file" << std::endl;
      ret = false;
    }
    else if (!checkOutputFileValid(opt.output))
    {
      std::cerr << "Error: unable to open output file" << std::endl;
      ret = false;
    }
  }
  
  if (opt.files.size() == 0)
  {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  }
  else
  {
    if (!opt.stream_in)
    {
      for (const auto &it : opt.files)
      {
        if (!checkFileExists(it))
        {
          std::cerr << "Error: File not found, " << it << std::endl;
          ret = false;
        }
      }
    }
  }
  
  return ret;
}

bool check_ProgramOptions_collapse(Bustools_opt& opt) {
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  // check for output directory
  if (opt.output.empty()) {
    std::cerr << "Error: Missing output directory" << std::endl;
    ret = false;
  }
  else {
    bool isDir = false;
    if (checkDirectoryExists(opt.output)) {
      isDir = true;
    }
    else {
      if (opt.output.at(opt.output.size() - 1) == '/') {
        if (my_mkdir(opt.output.c_str(), 0777) == -1) {
          std::cerr << "Error: could not create directory " << opt.output << std::endl;
          ret = false;
        }
        else {
          isDir = true;
        }
      }
    }
    
    if (isDir) {
      opt.output += "output";
    }
  }
  
  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  }
  else {
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
    ret = false;
  }
  else {
    if (!checkFileExists(opt.count_genes)) {
      std::cerr << "Error: File not found " << opt.count_genes << std::endl;
      ret = false;
    }
  }
  
  if (opt.count_ecs.size() == 0) {
    std::cerr << "Error: missing equivalence class mapping file" << std::endl;
    ret = false;
  }
  else {
    if (!checkFileExists(opt.count_ecs)) {
      std::cerr << "Error: File not found " << opt.count_ecs << std::endl;
      ret = false;
    }
  }
  
  if (opt.count_txp.size() == 0) {
    std::cerr << "Error: missing transcript name file" << std::endl;
    ret = false;
  }
  else {
    if (!checkFileExists(opt.count_txp)) {
      std::cerr << "Error: File not found " << opt.count_txp << std::endl;
      ret = false;
    }
  }
  
  return ret;
}

bool check_ProgramOptions_clusterhist(Bustools_opt& opt) {
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  // check for output directory
  if (opt.output.empty()) {
    std::cerr << "Error: Missing output directory" << std::endl;
    ret = false;
  }
  else {
    bool isDir = false;
    if (checkDirectoryExists(opt.output)) {
      isDir = true;
    }
    else {
      if (opt.output.at(opt.output.size() - 1) == '/') {
        if (my_mkdir(opt.output.c_str(), 0777) == -1) {
          std::cerr << "Error: could not create directory " << opt.output << std::endl;
          ret = false;
        }
        else {
          isDir = true;
        }
      }
    }
    
    std::string histDir = opt.output + "cluster_hists/";
    //generate directory
    if (!checkDirectoryExists(histDir)) {
      if (my_mkdir(histDir.c_str(), 0777) == -1) {
        std::cerr << "Error: could not create directory " << opt.output << std::endl;
        ret = false;
      }
    }
  }
  
  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing BUS input files" << std::endl;
    ret = false;
  }
  else {
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
    ret = false;
  }
  else {
    if (!checkFileExists(opt.count_genes)) {
      std::cerr << "Error: File not found " << opt.count_genes << std::endl;
      ret = false;
    }
  }
  
  if (opt.count_ecs.size() == 0) {
    std::cerr << "Error: missing equivalence class mapping file" << std::endl;
    ret = false;
  }
  else {
    if (!checkFileExists(opt.count_ecs)) {
      std::cerr << "Error: File not found " << opt.count_ecs << std::endl;
      ret = false;
    }
  }
  
  if (opt.count_txp.size() == 0) {
    std::cerr << "Error: missing transcript name file" << std::endl;
    ret = false;
  }
  else {
    if (!checkFileExists(opt.count_txp)) {
      std::cerr << "Error: File not found " << opt.count_txp << std::endl;
      ret = false;
    }
  }
  
  if (opt.cluster_input_file.size() == 0) {
    std::cerr << "Error: missing cluster file" << std::endl;
    ret = false;
  }
  else {
    if (!checkFileExists(opt.cluster_input_file)) {
      std::cerr << "Error: File not found " << opt.cluster_input_file << std::endl;
      ret = false;
    }
  }
  
  return ret;
}

bool check_ProgramOptions_extract(Bustools_opt &opt)
{
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  if (opt.output.empty())
  {
    std::cerr << "Error: missing output directory" << std::endl;
  }
  else
  {
    // check if output directory exists or if we can create it
    struct stat stFileInfo;
    auto intStat = stat(opt.output.c_str(), &stFileInfo);
    if (intStat == 0)
    {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode))
      {
        std::cerr << "Error: file " << opt.output << " exists and is not a directory" << std::endl;
        ret = false;
      }
    }
    else
    {
      // create directory
      if (my_mkdir(opt.output.c_str(), 0777) == -1)
      {
        std::cerr << "Error: could not create directory " << opt.output << std::endl;
        ret = false;
      }
    }
  }
  
  if (opt.files.size() == 0)
  {
    std::cerr << "Error: Missing BUS input file" << std::endl;
    ret = false;
  }
  else if (opt.files.size() == 1)
  {
    if (!opt.stream_in)
    {
      for (const auto &it : opt.files)
      {
        if (!checkFileExists(it))
        {
          std::cerr << "Error: File not found, " << it << std::endl;
          ret = false;
        }
      }
    }
  }
  else
  {
    std::cerr << "Error: Only one input file allowed" << std::endl;
    ret = false;
  }
  
  if (opt.fastq.size() == 0)
  {
    std::cerr << "Error: Missing FASTQ file(s)" << std::endl;
    ret = false;
  }
  else
  {
    for (const auto &f : opt.fastq)
    {
      if (!checkFileExists(f))
      {
        std::cerr << "Error: File not found, " << f << std::endl;
        ret = false;
      }
    }
  }
  
  if (opt.nFastqs == 0)
  {
    std::cerr << "Error: nFastqs is zero" << std::endl;
    ret = false;
  }
  else
  {
    if (opt.fastq.size() % opt.nFastqs != 0)
    {
      std::cerr << "Error: incorrect number of FASTQ file(s)" << std::endl;
      ret = false;
    }
  }

  if (opt.extract_exclude && opt.extract_include)
  {
    std::cerr << "Error: cannot specify both --exclude and --include" << std::endl;
    ret = false;
  }
  
  return ret;
}

bool check_ProgramOptions_inflate(Bustools_opt &opt)
{
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  if (!opt.stream_out)
  {
    if (opt.output.empty())
    {
      std::cerr << "Error: missing output file" << std::endl;
      ret = false;
    }
    else if (!checkOutputFileValid(opt.output))
    {
      std::cerr << "Error: unable to open output file" << std::endl;
      ret = false;
    }
  }

  if (opt.files.size() != 1)
  {
    ret = false;
    if (opt.files.size() == 0)
    {
      std::cerr << "Error: Missing BUSZ input file" << std::endl;
    }
    else
    {
      std::cerr << "Error: Multiple files not yet supported" << std::endl;
    }
  }
  else if (!opt.stream_in)
  {
    if (!checkFileExists(opt.files[0]))
    {
      std::cerr << "Error: File not found, " << opt.files[0] << std::endl;
      ret = false;
    }
  }

  return ret;
}

bool check_ProgramOptions_compress(Bustools_opt &opt)
{
  bool ret = true;

  if (opt.parse_error) {
    ret = false;
  }

  if (!opt.stream_out)
  {
    if (opt.output.empty())
    {
      std::cerr << "Error: missing output file" << std::endl;
      ret = false;
    }
    else if (!checkOutputFileValid(opt.output))
    {
      std::cerr << "Error: unable to open output file" << std::endl;
      ret = false;
    }
  }

  if (opt.files.size() != 1)
  {
    ret = false;
    if (opt.files.size() == 0)
    {
      std::cerr << "Error: Missing BUS input file" << std::endl;
    }
    else
    {
      std::cerr << "Error: Multiple files not yet supported" << std::endl;
    }
  }
  else if (!opt.stream_in)
  {
    if (!checkFileExists(opt.files[0]))
    {
      std::cerr << "Error: File not found, " << opt.files[0] << std::endl;
      ret = false;
    }
  }

  return ret;
}

void Bustools_Usage()
{
  std::cout << "bustools " << BUSTOOLS_VERSION << std::endl
            << std::endl
            << "Usage: bustools <CMD> [arguments] .." << std::endl
            << std::endl
            << "Where <CMD> can be one of: " << std::endl
            << std::endl
            << "sort            Sort a BUS file by barcodes and UMIs" << std::endl
            << "correct         Error correct a BUS file" << std::endl
            << "umicorrect      Error correct the UMIs in a BUS file" << std::endl
            << "count           Generate count matrices from a BUS file" << std::endl
            << "inspect         Produce a report summarizing a BUS file" << std::endl
            << "allowlist       Generate an on-list from a BUS file" << std::endl
            //<< "project         Project a BUS file to gene sets" << std::endl
            << "capture         Capture records from a BUS file" << std::endl
            //<< "merge           Merge bus files from same experiment" << std::endl
            << "text            Convert a binary BUS file to a tab-delimited text file" << std::endl
            << "fromtext        Convert a tab-delimited text file to a binary BUS file" << std::endl
            << "extract         Extract FASTQ reads correspnding to reads in BUS file" << std::endl
            //<< "predict         Correct the count matrix using prediction of unseen species" << std::endl
            //<< "collapse        Turn BUS files into a BUG file" << std::endl
            //<< "clusterhist     Create UMI histograms per cluster" << std::endl
            //<< "linker          Remove section of barcodes in BUS files" << std::endl
            << "compress        Compress a BUS file" << std::endl
            << "decompress      Decompress a BUSZ (compressed BUS) file" << std::endl
            << "version         Prints version number" << std::endl
            << "cite            Prints citation information" << std::endl
            << std::endl
            << "Running bustools <CMD> without arguments prints usage information for <CMD>"
            << std::endl
            << std::endl;
}

void Bustools_sort_Usage()
{
  std::cout << "Usage: bustools sort [options] bus-files" << std::endl
            << std::endl
            << "Options: " << std::endl
            << "Default behavior is to sort by barcode, UMI, ec, then flag" << std::endl
            << "-t, --threads         Number of threads to use" << std::endl
            << "-m, --memory          Maximum memory used" << std::endl
            << "-T, --temp            Location and prefix for temporary files " << std::endl
            << "                      required if using -p, otherwise defaults to output" << std::endl
            << "-o, --output          File for sorted output" << std::endl
            << "-p, --pipe            Write to standard output" << std::endl
            << "    --umi             Sort by UMI, barcode, then ec" << std::endl
            << "    --count           Sort by multiplicity, barcode, UMI, then ec" << std::endl
            << "    --flags           Sort by flag, ec, barcode, then UMI" << std::endl
            << "    --flags-bc        Sort by flag, barcode, UMI, then ec" << std::endl
            << "    --no-flags        Ignore and reset the flag while sorting" << std::endl
            << std::endl;
}

void Bustools_capture_Usage()
{
  std::cout << "Usage: bustools capture [options] bus-files" << std::endl
            << std::endl
            << "Options: " << std::endl
            << "-o, --output          File for captured output " << std::endl
            << "-x, --complement      Take complement of captured set" << std::endl
            << "-c, --capture         Capture list" << std::endl
            << "-e, --ecmap           File for mapping equivalence classes to transcripts (for --transcripts)" << std::endl
            << "-t, --txnames         File with names of transcripts (for --transcripts)" << std::endl
            << "-p, --pipe            Write to standard output" << std::endl
            << "Capture types: " << std::endl
            << "-F, --flags           Capture list is a list of flags to capture" << std::endl
            << "-s, --transcripts     Capture list is a list of transcripts to capture" << std::endl
            << "-u, --umis            Capture list is a list of UMIs to capture" << std::endl
            << "-b, --barcode         Capture list is a list of barcodes to capture" << std::endl
            << std::endl;
}

void Bustools_merge_Usage()
{
  std::cout << "Usage: bustools merge [options] sorted-bus-file by flag" << std::endl
            << std::endl
            << "Options: " << std::endl
            << "-o, --output          Directory for merged output" << std::endl
            << "-e, --ecmap           File for mapping equivalence classes to transcripts" << std::endl
            << "-t, --txnames         File with names of transcripts" << std::endl
            << std::endl;
}

void Bustools_mash_Usage()
{
  std::cout << "Usage: bustools mash [options] directories" << std::endl
            << "  Note: BUS files should be sorted by flag" << std::endl
            << std::endl
            << "Options: " << std::endl
            << "-t, --threads         Number of threads to use" << std::endl
            << "-o, --output          Directory for mashed output" << std::endl
            << std::endl;
}

void Bustools_dump_Usage()
{
  std::cout << "Usage: bustools text [options] bus-files" << std::endl
            << std::endl
            << "Options: " << std::endl
            << "-o, --output          File for text output" << std::endl
            << "-f, --flags           Write the flag column" << std::endl
            << "-d, --pad             Write the pad column" << std::endl
            << "-p, --pipe            Write to standard output" << std::endl
            << "-a, --showAll         Show hidden metadata in barcodes" << std::endl
            << std::endl;
}

void Bustools_fromtext_Usage() 
{
  std::cout << "Usage: bustools fromtext [options] text-files" << std::endl << std::endl
            << "Options: " << std::endl
            << "-o, --output          File for BUS output" << std::endl
            << "-p, --pipe            Write to standard output" << std::endl
            << std::endl;
}

void Bustools_correct_Usage()
{
  std::cout << "Usage: bustools correct [options] bus-files" << std::endl
            << std::endl
            << "Options: " << std::endl
            << "-o, --output          File for corrected bus output" << std::endl
            << "-w, --onlist          File of on-list barcodes to correct to" << std::endl
            << "-p, --pipe            Write to standard output" << std::endl
            << "-d, --dump            Dump uncorrected to corrected barcodes (optional)" << std::endl
            << "-r, --replace         The file of on-list barcodes is a barcode replacement file" << std::endl
            << "    --nocorrect       Skip barcode error correction and only keep perfect matches to on-list" << std::endl
            << std::endl;
}

void Bustools_count_Usage()
{
  std::cout << "Usage: bustools count [options] sorted-bus-files" << std::endl
            << std::endl
            << "Options: " << std::endl
            << "-o, --output          Output directory gene matrix files" << std::endl
            << "-g, --genemap         File for mapping transcripts to genes" << std::endl
            << "-e, --ecmap           File for mapping equivalence classes to transcripts" << std::endl
            << "-t, --txnames         File with names of transcripts" << std::endl
            << "    --genecounts      Aggregate counts to genes only" << std::endl
            << "    --umi-gene        Perform gene-level collapsing of UMIs" << std::endl
            //<< "    --em              Estimate gene abundances using EM algorithm" << std::endl
            << "    --cm              Count multiplicities instead of UMIs" << std::endl
            << "-s, --split           Split output matrix in two (plus ambiguous) based on transcripts supplied in this file" << std::endl
            << "-m, --multimapping    Include bus records that pseudoalign to multiple genes" << std::endl
            //<< "    --hist            Output copy per UMI histograms for all genes" << std::endl 
            //<< "-d  --downsample      Specify a factor between 0 and 1 specifying how much to downsample" << std::endl 
            //<< "    --rawcounts       The count matrix will contain raw counts instead of UMI counts" << std::endl 
            << std::endl;
}

void Bustools_predict_Usage() {
  std::cout << "Usage: bustools predict [options] count_output_dir" << std::endl << std::endl
            << "Options: " << std::endl
            << "-o, --output          Output directory" << std::endl
            << "-t, --predict_t       Specifies prediction range, will at predict t times the number of reads" << std::endl
            << std::endl;
}

void Bustools_umicorrect_Usage() 
{
  std::cout << "Usage: bustools umicorrect [options] sorted-bus-files" << std::endl << std::endl
            << "Options: " << std::endl
            << "-o, --output          Output directory gene matrix files" << std::endl
            << "-p, --pipe            Write to standard output" << std::endl
            << "-g, --genemap         File for mapping transcripts to genes" << std::endl
            << "-e, --ecmap           File for mapping equivalence classes to transcripts" << std::endl
            << "-t, --txnames         File with names of transcripts" << std::endl
            << std::endl;
}

void Bustools_whitelist_Usage()
{
  std::cout << "Usage: bustools allowlist [options] sorted-bus-file" << std::endl
            << std::endl
            << "Options: " << std::endl
            << "-o, --output        File for the on-list" << std::endl
            << "-f, --threshold     Minimum number of times a barcode must appear to be included in on-list" << std::endl
            << std::endl;
}

void Bustools_project_Usage()
{
  std::cout << "Usage: bustools project [options] sorted-bus-file" << std::endl
            << std::endl
            << "Options: " << std::endl
            << "-o, --output          File for project bug output and list of genes (no extension)" << std::endl
            << "-m, --map             File for mapping source to destination" << std::endl
            << "-e, --ecmap           File for mapping equivalence classes to transcripts" << std::endl
            << "-t, --txnames         File with names of transcripts" << std::endl
            << "-p, --pipe            Write to standard output" << std::endl
            << "Project types: " << std::endl
            << "-F, --flags           Map is a two column list of source flags to destination flags " << std::endl
            << "-s, --transcripts     [Output folder needed] Map is a two column list of source transcripts to destination transcripts" << std::endl
            << "-u, --umis            Map is a two column list of source umis to destination umis" << std::endl
            << "-b, --barcode         Map is a two column list of source barcodes to destination barcodes" << std::endl
            << std::endl;
}

void Bustools_inspect_Usage()
{
  std::cout << "Usage: bustools inspect [options] sorted-bus-file" << std::endl
            << std::endl
            << "Options: " << std::endl
            << "-o, --output          File for JSON output (optional)" << std::endl
            << "-e, --ecmap           File for mapping equivalence classes to transcripts" << std::endl
            << "-w, --onlist          File of on-list barcodes to correct to" << std::endl
            << "-p, --pipe            Write to standard output" << std::endl
            << std::endl;
}

void Bustools_linker_Usage()
{
  std::cout << "Usage: bustools linker [options] bus-files" << std::endl
            << std::endl
            << "Options: " << std::endl
            << "-o, --output          Output BUS file" << std::endl
            << "-s, --start           Start coordinate for section of barcode to remove (0-indexed, inclusive)" << std::endl
            << "-e, --end             End coordinate for section of barcode to remove (0-indexed, exclusive)" << std::endl
            << "-p, --pipe            Write to standard output" << std::endl
            << std::endl;
}

void Bustools_collapse_Usage() 
{
  std::cout << "Usage: bustools collapse [options] sorted-bus-files" << std::endl << std::endl
            << "Options: " << std::endl
            << "-o, --output          Output directory gene matrix files" << std::endl
            << "-g, --genemap         File for mapping transcripts to genes" << std::endl
            << "-e, --ecmap           File for mapping equivalence classes to transcripts" << std::endl
            << "-t, --txnames         File with names of transcripts" << std::endl
            << "-p, --pipe            Write to standard output" << std::endl
            << std::endl;
}

void Bustools_clusterhist_Usage() 
{
  std::cout << "Usage: bustools clusterhist [options] sorted-bus-files" << std::endl << std::endl
            << "Options: " << std::endl
            << "-o, --output          Output directory gene matrix files" << std::endl
            << "-g, --genemap         File for mapping transcripts to genes" << std::endl
            << "-e, --ecmap           File for mapping equivalence classes to transcripts" << std::endl
            << "-t, --txnames         File with names of transcripts" << std::endl
            << "-c, --clusterfile     File with cell cluster assignments" << std::endl
            << "-p, --pipe            Write to standard output" << std::endl
            << std::endl;
}

void Bustools_extract_Usage()
{
  std::cout << "Usage: bustools extract [options] sorted-bus-file" << std::endl
            << "  Note: BUS file should be sorted by flag using bustools sort --flag" << std::endl
            << std::endl
            << "Options: " << std::endl
            << "-o, --output          Output directory for FASTQ files" << std::endl
            << "-f, --fastq           FASTQ file(s) from which to extract reads (comma-separated list)" << std::endl
            << "-N, --nFastqs         Number of FASTQ file(s) per run" << std::endl
            << "-x, --exclude         Exclude reads in the BUS file from the specified FASTQ file(s)" << std::endl
            << "-i, --include         Include reads in the BUS file from the specified FASTQ file(s)" << std::endl
            << std::endl;
}

void Bustools_compress_Usage()
{
  std::cout << "Usage: bustools compress [options] sorted-bus-file" << std::endl
            << "Note: BUS file should be sorted by barcode-umi-ec" << std::endl
            << std::endl
            << "Options: " << std::endl
            << "-N, --chunk-size CHUNK_SIZE    Number of rows to compress as a single block." << std::endl
            // << "-L, --lossy-umi                Allow lossy compression over UMIs. Each UMI will be renamed for minimal compression." << std::endl
            << "-o, --output OUTPUT            Write compressed file to OUTPUT." << std::endl
            << "-p, --pipe                     Write to standard output." << std::endl
            << "-h, --help                     Print this message and exit." << std::endl
            << std::endl;
}

void Bustools_inflate_Usage()
{
  std::cout << "Usage: bustools {inflate | decompress} [options] compressed-bus-file" << std::endl
            << std::endl
            << "Options: " << std::endl
            << "-p, --pipe               Write to standard output." << std::endl
            << "-o, --output OUTPUT      File for inflated output." << std::endl
            << "-h, --help               Print this message and exit." << std::endl
            << std::endl;
}

void print_citation()
{
  std::cout << "When using this program in your research, please cite" << std::endl
            << std::endl
            << "  Melsted, P., Booeshaghi, A. S., et al." << std::endl
            << "  Modular, efficient and constant-memory single-cell RNA-seq preprocessing, " << std::endl
            << "  Nature Biotechnology (2021), doi:10.1038/s41587-021-00870-2" << std::endl
            << std::endl;
}

void print_version()
{
  std::cout << "bustools, version " << BUSTOOLS_VERSION << std::endl;
}

int main(int argc, char **argv)
{
  
  if (argc < 2)
  {
    // Print error message, function?
    Bustools_Usage();
    exit(1);
  }
  else
  {
    bool disp_help = argc == 2;
    std::string cmd(argv[1]);
    Bustools_opt opt;
    
    if (cmd == "sort")
    {
      if (disp_help)
      {
        Bustools_sort_Usage();
        exit(0);
      }
      parse_ProgramOptions_sort(argc - 1, argv + 1, opt);
      if (check_ProgramOptions_sort(opt))
      { //Program options are valid
        //bustools_sort_orig(opt);
        bustools_sort(opt);
      }
      else
      {
        Bustools_sort_Usage();
        exit(1);
      }
    }
    else if (cmd == "merge")
    {
      if (disp_help)
      {
        Bustools_merge_Usage();
        exit(0);
      }
      parse_ProgramOptions_merge(argc - 1, argv + 1, opt);
      if (check_ProgramOptions_merge(opt))
      {
        bustools_merge_different_index(opt);
        // bustools_merge(opt);
      }
      else
      {
        Bustools_merge_Usage();
        exit(1);
      }
    }
    else if (cmd == "mash")
    {
      if (disp_help)
      {
        Bustools_mash_Usage();
        exit(0);
      }
      parse_ProgramOptions_mash(argc - 1, argv + 1, opt);
      if (check_ProgramOptions_mash(opt))
      {
        bustools_mash(opt);
        // bustools_mash(opt);
      }
      else
      {
        Bustools_mash_Usage();
        exit(1);
      }
    }
    else if (cmd == "cite")
    {
      print_citation();
    }
    else if (cmd == "version")
    {
      print_version();
    }
    else if (cmd == "dump" || cmd == "text")
    {
      if (disp_help)
      {
        Bustools_dump_Usage();
        exit(0);
      }
      parse_ProgramOptions_dump(argc - 1, argv + 1, opt);
      if (check_ProgramOptions_dump(opt))
      { //Program options are valid
        bustools_text(opt);
      }
      else
      {
        Bustools_dump_Usage();
        exit(1);
      }
    }
    else if (cmd == "fromtext") {
      if (disp_help) 
      {
        Bustools_fromtext_Usage();
        exit(0);        
      }
      parse_ProgramOptions_fromtext(argc-1, argv+1, opt);
      if (check_ProgramOptions_fromtext(opt)) 
      { 
        //Program options are valid
        bustools_fromtext(opt); //found in the bustools_text.cpp file
      } 
      else 
      {
        Bustools_fromtext_Usage();
        exit(1);
      }
    }
    else if (cmd == "correct")
    {
      if (disp_help)
      {
        Bustools_correct_Usage();
        exit(0);
      }
      parse_ProgramOptions_correct(argc - 1, argv + 1, opt);
      if (check_ProgramOptions_correct(opt))
      { //Program options are valid
        if (opt.split_correct)
        {
          bustools_split_correct(opt);
        }
        else
        {
          bustools_correct(opt);
        }
      }
      else
      {
        Bustools_dump_Usage();
        exit(1);
      }
    }
    else if (cmd == "count")
    {
      if (disp_help)
      {
        Bustools_count_Usage();
        exit(0);
      }
      parse_ProgramOptions_count(argc - 1, argv + 1, opt);
      if (check_ProgramOptions_count(opt))
      { //Program options are valid
        if (opt.count_cm)
        {
          bustools_count_mult(opt);
        }
        else
        {
          bustools_count(opt);
        }
      }
      else
      {
        Bustools_count_Usage();
        exit(1);
      }
    }
    else if (cmd == "predict") 
    {
      if (disp_help) 
      {
        Bustools_predict_Usage();
        exit(0);        
      }
      parse_ProgramOptions_predict(argc-1, argv+1, opt);
      if (check_ProgramOptions_predict(opt)) 
      { 
        bustools_predict(opt);
      } else {
        Bustools_predict_Usage();
        exit(1);
      }
    } 
    else if (cmd == "umicorrect") 
    {
      if (disp_help) 
      {
        Bustools_umicorrect_Usage();
        exit(0);        
      }
      parse_ProgramOptions_umicorrect(argc-1, argv+1, opt);
      if (check_ProgramOptions_umicorrect(opt)) 
      {
        bustools_umicorrect(opt);
      } 
      else 
      {
        Bustools_umicorrect_Usage();
        exit(1);
      }
    }
    else if (cmd == "capture")
    {
      if (disp_help)
      {
        Bustools_capture_Usage();
        exit(0);
      }
      parse_ProgramOptions_capture(argc - 1, argv + 1, opt);
      if (check_ProgramOptions_capture(opt))
      { //Program options are valid
        bustools_capture(opt);
      }
      else
      {
        Bustools_capture_Usage();
        exit(1);
      }
    }
    else if (cmd == "whitelist" || cmd == "allowlist" || cmd == "onlist")
    {
      if (disp_help)
      {
        Bustools_whitelist_Usage();
        exit(0);
      }
      parse_ProgramOptions_whitelist(argc - 1, argv + 1, opt);
      if (check_ProgramOptions_whitelist(opt))
      { //Program options are valid
        bustools_whitelist(opt);
      }
      else
      {
        Bustools_whitelist_Usage();
        exit(1);
      }
    }
    else if (cmd == "project")
    {
      if (disp_help)
      {
        Bustools_project_Usage();
        exit(0);
      }
      parse_ProgramOptions_project(argc - 1, argv + 1, opt);
      if (check_ProgramOptions_project(opt))
      { //Program options are valid
        bustools_project(opt);
      }
      else
      {
        Bustools_project_Usage();
        exit(1);
      }
    }
    else if (cmd == "inspect")
    {
      if (disp_help)
      {
        Bustools_inspect_Usage();
        exit(0);
      }
      parse_ProgramOptions_inspect(argc - 1, argv + 1, opt);
      if (check_ProgramOptions_inspect(opt))
      { //Program options are valid
        bustools_inspect(opt);
      }
      else
      {
        Bustools_inspect_Usage();
        exit(1);
      }
    }
    else if (cmd == "linker")
    {
      if (disp_help)
      {
        Bustools_linker_Usage();
        exit(0);
      }
      parse_ProgramOptions_linker(argc - 1, argv + 1, opt);
      if (check_ProgramOptions_linker(opt))
      { //Program options are valid
        bustools_linker(opt);
      }
      else
      {
        Bustools_linker_Usage();
        exit(1);
      }
    }
    else if (cmd == "collapse") 
    {
      if (disp_help) 
      {
        Bustools_collapse_Usage();
        exit(0);        
      }
      parse_ProgramOptions_collapse(argc-1, argv+1, opt);
      if (check_ProgramOptions_collapse(opt)) 
      {
        bustools_collapse(opt);
      } 
      else 
      {
        Bustools_collapse_Usage();
        exit(1);
      }
    } 
    else if (cmd == "clusterhist") 
    {
      if (disp_help) {
        Bustools_clusterhist_Usage();
        exit(0);        
      }
      parse_ProgramOptions_clusterhist(argc-1, argv+1, opt);
      if (check_ProgramOptions_clusterhist(opt)) 
      {
        bustools_clusterhist(opt);
      } 
      else 
      {
        Bustools_clusterhist_Usage();
        exit(1);
      }
    }
    else if (cmd == "extract")
    {
      if (disp_help)
      {
        Bustools_extract_Usage();
        exit(0);
      }
      parse_ProgramOptions_extract(argc - 1, argv + 1, opt);
      if (check_ProgramOptions_extract(opt))
      { //Program options are valid
        bustools_extract(opt);
      }
      else
      {
        Bustools_extract_Usage();
        exit(1);
      }
    }
    else if(cmd == "compress")
    {
      if (disp_help)
      {
        Bustools_compress_Usage();
        exit(0);
      }
      if (parse_ProgramOptions_compress(argc - 1, argv + 1, opt))
      {
        Bustools_compress_Usage();
        exit(0);
      }
      else if (check_ProgramOptions_compress(opt))
      {
        bustools_compress(opt);
        exit(0);
      }
      else
      {
        Bustools_compress_Usage();
        exit(1);
      }
    }
    else if(cmd == "inflate" || cmd == "decompress")
    {
      if (disp_help)
      {
        Bustools_inflate_Usage();
        exit(0);
      }
      if (parse_ProgramOptions_inflate(argc - 1, argv + 1, opt))
      {
        Bustools_inflate_Usage();
        exit(0);
      }
      else if (check_ProgramOptions_inflate(opt))
      {
        bustools_decompress(opt);
        exit(0);
      }
      else
      {
        Bustools_inflate_Usage();
        exit(1);
      }
    }
    else
    {
      std::cerr << "Error: invalid command " << cmd << std::endl;
      Bustools_Usage();
    }
  }
}
