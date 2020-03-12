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
#include "bustools_text.h"
#include "bustools_collapse.h"
#include "bustools_umicorrect.h"
#include "bustools_predquant.h"


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

bool checkOutputFileValid(const std::string &fn) {
  std::ofstream of(fn);
  if (of.is_open()) {
    of.close();
    return checkFileExists(fn);
  } else {
    return false;
  }
}

std::vector<std::string> parseList(const std::string &s, const std::string &sep = ",") {
  std::vector<std::string> ret;
  size_t start = 0, end;
  while ((end = s.find(sep, start)) != std::string::npos) {
    ret.push_back(s.substr(start, end - start));
    start = end + sep.length();
  }
  end = s.size();
  if (end - start > 0) {
    ret.push_back(s.substr(start, end - start));
  }
  return ret;
}



void parse_ProgramOptions_sort(int argc, char **argv, Bustools_opt& opt) {

  const char* opt_string = "t:o:m:T:cusp";

  static struct option long_options[] = {
    {"threads",         required_argument,  0, 't'},
    {"output",          required_argument,  0, 'o'},
    {"memory",          required_argument,  0, 'm'},
    {"temp",            required_argument,  0, 'T'},
    {"umi",             no_argument,        0, 'u'},
    {"count",           no_argument,        0, 'c'},
    {"flags",           no_argument,        0, 'F'},
    {"pipe",            no_argument,        0, 'p'},
    {0,                 0,                  0,  0 }
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {
    std::string s;
    size_t sh = 0;
    int n = 0;
    switch (c) {

    case 't':
      opt.threads = atoi(optarg);
      break;    
    case 'o':
      opt.output = optarg;
      break;
    case 'm':
      s = optarg;
      sh = 0;
      n = s.size();
      if (n==0) {
        break;
      }
      switch(s[n-1]) {
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
      opt.max_memory = atoi(s.substr(0,n).c_str());
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

void parse_ProgramOptions_capture(int argc, char **argv, Bustools_opt& opt) {
   const char* opt_string = "o:xc:e:t:Fsubfp";

  static struct option long_options[] = {
    {"output",          required_argument,  0, 'o'},
    {"complement",      no_argument,        0, 'x'},
    {"capture",         required_argument,  0, 'c'},
    {"ecmap",           required_argument,  0, 'e'},
    {"txnames",         required_argument,  0, 't'},
    {"flags",           no_argument,        0, 'F'},
    {"transcripts",     no_argument,        0, 's'},
    {"umis",            no_argument,        0, 'u'},
    {"barcode",         no_argument,        0, 'b'},
    {"combo",           no_argument,        0, 'f'},
    {"pipe",            no_argument,        0, 'p'},
    {0,                 0,                  0,  0 }
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

    switch (c) {
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
    default:
      break;
    }
  }

  while (optind < argc) opt.files.push_back(argv[optind++]);

  if (opt.files.size() == 1 && opt.files[0] == "-") {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_count(int argc, char **argv, Bustools_opt& opt) {
  const char* opt_string = "o:g:e:t:m";
  int gene_flag = 0;
  int em_flag = 0;
  static struct option long_options[] = {
    {"output",          required_argument,  0, 'o'},
    {"genemap",          required_argument,  0, 'g'},
    {"ecmap",          required_argument,  0, 'e'},
    {"txnames",          required_argument,  0, 't'},
    {"genecounts", no_argument, &gene_flag, 1},
    {"multimapping", no_argument, 0, 'm'},
    {"em", no_argument, &em_flag, 1},
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
    case 'm':
      opt.count_gene_multimapping = true;
      break;
    default:
      break;
    }
  }
  if (gene_flag) {
    opt.count_collapse = true;
  }
  if (em_flag) {
    opt.count_em = true;
  }

  while (optind < argc) opt.files.push_back(argv[optind++]);

  if (opt.files.size() == 1 && opt.files[0] == "-") {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_dump(int argc, char **argv, Bustools_opt& opt) {

  const char* opt_string = "o:pf";

  static struct option long_options[] = {
    {"output",          required_argument,  0, 'o'},
    {"pipe",            no_argument, 0, 'p'},
    {"flags",           no_argument, 0, 'f'},
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
    case 'f':
      opt.text_dumpflags = true;
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

void parse_ProgramOptions_fromtext(int argc, char **argv, Bustools_opt& opt) {

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

void parse_ProgramOptions_whitelist(int argc, char **argv, Bustools_opt &opt) {
  
  /* Parse options. */
  const char *opt_string = "o:f:h";

  static struct option long_options[] = {
    {"output", required_argument, 0, 'o'},
    {"threshold", required_argument, 0, 'f'},
    {0, 0, 0, 0}
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {
    switch (c) {
      case 'o':
        opt.output = optarg;
        break;
      case 'f':
        opt.threshold = atoi(optarg);
        break;
      default:
        break;
    }
  }

  /* All other argumuments are (sorted) BUS files. */
  while (optind < argc) opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-") {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_project(int argc, char **argv, Bustools_opt &opt) {
  
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
    {0, 0, 0, 0}
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {
    switch (c) {
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
      default:
        break;
    }
  }

  /* All other argumuments are (sorted) BUS files. */
  while (optind < argc) opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-") {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_inspect(int argc, char **argv, Bustools_opt &opt) {
  
  /* Parse options. */
  const char *opt_string = "o:e:w:p";

  static struct option long_options[] = {
    {"output", required_argument, 0, 'o'},
    {"ecmap", required_argument, 0, 'e'},
    {"whitelist", required_argument, 0, 'w'},
    {"pipe", no_argument, 0, 'p'},
    {0, 0, 0, 0}
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {
    switch (c) {
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
      default:
        break;
    }
  }

  /* All other argumuments are (sorted) BUS files. */
  while (optind < argc) opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-") {
    opt.stream_in = true;
  }
}

void parse_ProgramOptions_linker(int argc, char **argv, Bustools_opt &opt) {
  
  /* Parse options. */
  const char *opt_string = "o:s:e:p";

  static struct option long_options[] = {
    {"output", required_argument, 0, 'o'},
    {"start", required_argument, 0, 's'},
    {"end", required_argument, 0, 'e'},
    {"pipe", no_argument, 0, 'p'},
    {0, 0, 0, 0}
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {
    switch (c) {
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
      default:
        break;
    }
  }

  /* All other argumuments are (sorted) BUS files. */
  while (optind < argc) opt.files.push_back(argv[optind++]);
  
  if (opt.files.size() == 1 && opt.files[0] == "-") {
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
		default:
			break;
		}
	}

	while (optind < argc) opt.files.push_back(argv[optind++]);

	if (opt.files.size() == 1 && opt.files[0] == "-") {
		opt.stream_in = true;
	}
}

void parse_ProgramOptions_umicorrect(int argc, char** argv, Bustools_opt& opt) {

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

void parse_ProgramOptions_predquant(int argc, char** argv, Bustools_opt& opt) {
	const char* opt_string = "o:g:mGP:I:U:N:";
	int gt_flag = 0;
	static struct option long_options[] = {
	  {"output",          required_argument,  0, 'o'},
	  {"genenames",       required_argument,  0, 'g'},
	  {"multimapping",    no_argument,  0, 'm'},
	  {"goodtoulmin",     no_argument,  &gt_flag, 'G'},
	  {"predtarget",      required_argument, 0, 'P'},
	  {"incl_bucket_limit", required_argument, 0, 'I'},
	  {"use_bucket_limit", required_argument, 0, 'U'},
	  {"num_buckets", required_argument, 0, 'N'},
//	  {"em", no_argument, &em_flag, 1},
	  {0,                 0,                  0,  0 }
	};

	int option_index = 0, c;

	while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

		switch (c) {
		case 'o':
			opt.output = optarg;
			break;
		case 'g':
			opt.predquant_genefile = optarg;
			break;
		case 'm':
			opt.count_gene_multimapping = true;
			break;
		case 'P':
			opt.predquant_pred_target = atoi(optarg);
			break;
		case 'I':
			opt.predquant_include_bucket_limit = atoi(optarg);
			break;
		case 'U':
			opt.predquant_use_bucket_limit = atoi(optarg);
			break;
		case 'N':
			opt.predquant_num_buckets = atoi(optarg);
			break;
		default:
			break;
		}
	}
	if (gt_flag) {
		opt.predquant_algorithm = PREDQUANT_ALG::GOOD_TOULMIN;
	}
//	if (em_flag) {
//		opt.count_em = true;
//	}

	while (optind < argc) opt.files.push_back(argv[optind++]);

	if (opt.files.size() == 1 && opt.files[0] == "-") {
		opt.stream_in = true;
	}
}

void parse_ProgramOptions_extract(int argc, char **argv, Bustools_opt &opt) {
  
  /* Parse options. */
  const char *opt_string = "o:f:N:p";

  static struct option long_options[] = {
    {"output", required_argument, 0, 'o'},
    {"fastq", required_argument, 0, 'f'},
    {"nFastqs", required_argument, 0, 'N'},
    {"pipe", no_argument, 0, 'p'},
    {0, 0, 0, 0}
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {
    switch (c) {
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
      default:
        break;
    }
  }

  /* All other argumuments are (sorted) BUS files. */
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

  if (!opt.stream_out) {
    if (opt.output.empty()) {
      std::cerr << "Error: missing output file" << std::endl;
      ret = false;
    } else if (!checkOutputFileValid(opt.output)) {
      std::cerr << "Error: unable to open output file" << std::endl;
      ret = false;
    }
  } 

  if (opt.max_memory < 1ULL<<26) {
    if (opt.max_memory < 128) {
      std::cerr << "Warning: low number supplied for maximum memory usage without M or G suffix\n  interpreting this as " << opt.max_memory << "Gb" << std::endl;
      opt.max_memory <<= 30;
    } else {
      std::cerr << "Warning: low number supplied for maximum memory, defaulting to 64Mb" << std::endl;
      opt.max_memory = 1ULL<<26; // 64Mb is absolute minimum
    }
  }

  if (opt.temp_files.empty()) {
    if (opt.stream_out) {
      std::cerr << "Error: temporary location -T must be set when using streaming output" << std::endl;
      ret = false;      
    } else {
      opt.temp_files = opt.output + ".";
    }
  } else {
    
    if (checkDirectoryExists(opt.temp_files)) {
      // if it is a directory, create random file prefix
      opt.temp_files += "/bus.sort." + std::to_string(getpid()) + ".";
    } else {
      if (opt.temp_files.at(opt.temp_files.size()-1) == '/') {
        // try creating the directory
        if (my_mkdir(opt.temp_files.c_str(),0777) == 0) {
          // 
          opt.temp_files += "/bus.sort." + std::to_string(getpid()) + ".";
        } else {
          std::cerr << "Error: directory " << opt.temp_files << " does not exist and could not be created. Check that the parent directory exists and you have write permissions." << std::endl;
          ret = false;
        }  
      } else {
        int n = opt.temp_files.size();
        if (opt.temp_files[n-1] != '.') {
          opt.temp_files += '.';
        }
        // assume the directory exist       
      }
    }
  }

  if (ret) {
    FILE *tmpf = fopen(opt.temp_files.c_str(), "a+");
    if (tmpf == nullptr) {
      char *dirn = dirname(strdup(opt.temp_files.c_str()));
      std::cerr << "Error: Could not create a temporary file in directory " << dirn << ". check that the directory exists and if you have permissions to write in it" << std::endl;
      free(dirn);
      ret = false;
    }
    remove(opt.temp_files.c_str());
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
    std::cerr << "Error: missing output directory" << std::endl;
  } else {
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

  }

  return ret;
}

bool check_ProgramOptions_dump(Bustools_opt& opt) {
  bool ret = true;

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

bool check_ProgramOptions_fromtext(Bustools_opt& opt) {
  bool ret = true;

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

bool check_ProgramOptions_capture(Bustools_opt& opt) {
  bool ret = true;

  if (opt.filter) {
    // check if output directory exists or if we can create it
    if (!checkDirectoryExists(opt.output)) {
      std::cerr << "Error: file " << opt.output << " exists and is not a directory" << std::endl;
      ret = false;      
    } else {
      // create directory
      if (my_mkdir(opt.output.c_str(), 0777) == -1) {
        std::cerr << "Error: could not create directory " << opt.output << std::endl;
        ret = false;
      }
    }
  } else if (!opt.stream_out) {
    if (opt.output.empty()) {
      std::cerr << "Error: missing output file" << std::endl;
      ret = false;
    } else if (!checkOutputFileValid(opt.output)) {
      std::cerr << "Error: unable to open output file" << std::endl;
      ret = false;
    }
  }
 
  if (opt.capture.empty()) {
    std::cerr << "Error: missing capture list" << std::endl;
    ret = false;
  } else {
    if (!checkFileExists(opt.capture)) {
      std::cerr << "Error: File not found, " << opt.capture << std::endl;
      ret = false;
    }
  }

  if (opt.type == CAPTURE_NONE) {
    std::cerr << "Error: capture list type must be specified (one of -s, -u, or -b)" << std::endl;
    ret = false;
  }

  if (opt.type == CAPTURE_TX) {
    if (opt.count_ecs.size() == 0) {
      std::cerr << "Error: missing equialence class mapping file" << std::endl;
    } else {
      if (!checkFileExists(opt.count_ecs)) {
        std::cerr << "Error: File not found " << opt.count_ecs << std::endl;
        ret = false;
      }
    }

    if (opt.count_txp.size() == 0) {
      std::cerr << "Error: missing transcript name file" << std::endl;
    } else {
      if (!checkFileExists(opt.count_txp)) {
        std::cerr << "Error: File not found " << opt.count_txp << std::endl;
        ret = false;
      }
    }
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

  if (opt.filter && (opt.complement || opt.type != CAPTURE_TX)) {
    std::cerr << "Warning: filter only meaningful without complement flag, and to"
      << " capture transcripts; no new ec file will be generated" << std::endl;
    opt.filter = false;
  }

  return ret;
}

bool check_ProgramOptions_correct(Bustools_opt& opt) {
  bool ret = true;

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

  if (opt.count_em && opt.count_gene_multimapping) {
    std::cerr << "Error: EM algorithm and counting multimapping reads are incompatible" << std::endl;
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

bool check_ProgramOptions_whitelist(Bustools_opt &opt) {
  bool ret = true;

  if (opt.output.empty()) {
    std::cerr << "Error: missing output file" << std::endl;
    ret = false;
  } else if (!checkOutputFileValid(opt.output)) {
    std::cerr << "Error: unable to open output file" << std::endl;
    ret = false;
  }

  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing BUS input file" << std::endl;
    ret = false;
  } else if (opt.files.size() == 1) {
    if (!opt.stream_in) {
      for (const auto& it : opt.files) {
        if (!checkFileExists(it)) {
          std::cerr << "Error: File not found, " << it << std::endl;
          ret = false;
        }
      }
    }
  } else {
    std::cerr << "Error: Only one input file allowed" << std::endl;
    ret = false;
  }

  if (opt.threshold < 0) { // threshold = 0 for no threshold
    std::cerr << "Error: Threshold cannot be less than or equal to 0" << std::endl;
    ret = false;
  }

  return ret;
}

bool check_ProgramOptions_project(Bustools_opt &opt) {
  bool ret = true;

  if (opt.output.empty()) {
    std::cerr << "Error: Missing output file" << std::endl;
  } 
  else if (!checkOutputFileValid(opt.output)) {
    std::cerr << "Error: unable to open output file" << std::endl;
    ret = false;
  }
  if (opt.type == PROJECT_TX){
    if (opt.output_folder.empty()) {
      std::cerr << "Error: Missing output folder" << std::endl;
    } else {
      // check if output directory exists or if we can create it
      struct stat stFileInfo;
      auto intStat = stat(opt.output_folder.c_str(), &stFileInfo);
      if (intStat == 0) {
        // file/dir exits
        if (!S_ISDIR(stFileInfo.st_mode)) {
          std::cerr << "Error: file " << opt.output_folder << " exists and is not a directory" << std::endl;
          ret = false;
        } 
      } else {
        // create directory
        if (my_mkdir(opt.output_folder.c_str(), 0777) == -1) {
          std::cerr << "Error: could not create directory " << opt.output_folder << std::endl;
          ret = false;
        }
      }
    }
  }

  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing BUS input file" << std::endl;
    ret = false;
  } else if (opt.files.size() == 1) {
    if (!opt.stream_in) {
      for (const auto& it : opt.files) {
        if (!checkFileExists(it)) {
          std::cerr << "Error: File not found, " << it << std::endl;
          ret = false;
        }
      }
    }
  } else {
    std::cerr << "Error: Only one input file allowed" << std::endl;
    ret = false;
  }

  if (opt.map.size() == 0) {
    std::cerr << "Error: missing mapping file" << std::endl;
  } else {
    if (!checkFileExists(opt.map)) {
      std::cerr << "Error: File not found " << opt.map << std::endl;
      ret = false;
    }
  }
  
  if (opt.count_ecs.size() == 0) {
    std::cerr << "Error: missing equivalence class mapping file" << std::endl;
  } else {
    if (!checkFileExists(opt.count_ecs)) {
      std::cerr << "Error: File not found " << opt.count_ecs << std::endl;
      ret = false;
    }
  }
  
  if (opt.count_txp.size() == 0) {
    std::cerr << "Error: missing transcript name file" << std::endl;
  } else {
    if (!checkFileExists(opt.map)) {
      std::cerr << "Error: File not found " << opt.count_txp << std::endl;
      ret = false;
    }
  }

  return ret;
}

bool check_ProgramOptions_inspect(Bustools_opt &opt) {
  bool ret = true;

  if (opt.output.size() && !checkOutputFileValid(opt.output)) {
    std::cerr << "Error: unable to open output file" << std::endl;
    ret = false;
  }
  
  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing BUS input file" << std::endl;
    ret = false;
  } else if (opt.files.size() == 1) {
    if (!opt.stream_in) {
      for (const auto& it : opt.files) {
        if (!checkFileExists(it)) {
          std::cerr << "Error: File not found, " << it << std::endl;
          ret = false;
        }
      }
    }
  } else {
    std::cerr << "Error: Only one input file allowed" << std::endl;
    ret = false;
  }
  
  if (opt.count_ecs.size()) {
    if (!checkFileExists(opt.count_ecs)) {
      std::cerr << "Error: File not found " << opt.count_ecs << std::endl;
      ret = false;
    }
  }

  if (opt.whitelist.size()) {
    if (!checkFileExists(opt.whitelist)) {
      std::cerr << "Error: File not found " << opt.whitelist << std::endl;
      ret = false;
    }
  }
  
  return ret;
}

bool check_ProgramOptions_linker(Bustools_opt &opt) {
  bool ret = true;
  
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
  
  return ret;
}

bool check_ProgramOptions_collapse(Bustools_opt& opt) {
	bool ret = true;

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

bool check_ProgramOptions_umicorrect(Bustools_opt& opt) {
	bool ret = true;

	if (!opt.stream_out) {
		if (opt.output.empty()) {
			std::cerr << "Error: missing output file" << std::endl;
			ret = false;
		}
		else if (!checkOutputFileValid(opt.output)) {
			std::cerr << "Error: unable to open output file" << std::endl;
			ret = false;
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

	return ret;
}


bool check_ProgramOptions_predquant(Bustools_opt& opt) {
	bool ret = true;

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

//unclear if EM should be in here
//	if (opt.count_em && opt.count_gene_multimapping) {
//		std::cerr << "Error: EM algorithm and counting multimapping reads are incompatible" << std::endl;
//		ret = false;
//	}



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

	if (opt.predquant_genefile.size() == 0) {
		std::cerr << "Error: missing genenames file" << std::endl;
		ret = false;
	}
	else {
		if (!checkFileExists(opt.predquant_genefile)) {
			std::cerr << "Error: File not found " << opt.predquant_genefile << std::endl;
			ret = false;
		}
	}
	
	//Good-Toulmin is incompatible with predtarg, always predicts to 2 - not sure how to solve that though, skip this for now
	//TODO: Don't do any checks on these params for now, fix later
	//{"goodtoulmin", no_argument, 0, 'G'},
	//{ "predtarg",        required_argument,  0, 'P' },
	//{ "incl_bucket_limit",  required_argument, 0, 'I' },
	//{ "use_bucket_limit",  required_argument, 0, 'U' },
	//{ "num_buckets",     required_argument, 0, 'N' },

	return ret;
}

bool check_ProgramOptions_extract(Bustools_opt &opt) {
  bool ret = true;
  
  if (opt.output.empty()) {
    std::cerr << "Error: missing output directory" << std::endl;
  } else {
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

  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing BUS input file" << std::endl;
    ret = false;
  } else if (opt.files.size() == 1) {
    if (!opt.stream_in) {
      for (const auto& it : opt.files) {
        if (!checkFileExists(it)) {
          std::cerr << "Error: File not found, " << it << std::endl;
          ret = false;
        }
      }
    }
  } else {
    std::cerr << "Error: Only one input file allowed" << std::endl;
    ret = false;
  }

  if (opt.fastq.size() == 0) {
    std::cerr << "Error: Missing FASTQ file(s)" << std::endl;
    ret = false;
  } else {
    for (const auto &f : opt.fastq) {
      if (!checkFileExists(f)) {
        std::cerr << "Error: File not found, " << f << std::endl;
        ret = false;
      }
    }
  }

  if (opt.nFastqs == 0) {
    std::cerr << "Error: nFastqs is zero" << std::endl;
    ret = false;
  } else {
    if (opt.fastq.size() % opt.nFastqs != 0) {
      std::cerr << "Error: incorrect number of FASTQ file(s)" << std::endl;
      ret = false;
    }
  }
  
  return ret;
}


void Bustools_Usage() {
  std::cout << "bustools " << BUSTOOLS_VERSION << std::endl << std::endl  
  << "Usage: bustools <CMD> [arguments] .." << std::endl << std::endl
  << "Where <CMD> can be one of: " << std::endl << std::endl
  << "sort            Sort a BUS file by barcodes and UMIs" << std::endl
  << "correct         Error correct a BUS file" << std::endl
  << "count           Generate count matrices from a BUS file" << std::endl
  << "inspect         Produce a report summarizing a BUS file" << std::endl
  << "whitelist       Generate a whitelist from a BUS file" << std::endl
  << "project         Project a BUS file to gene sets" << std::endl
  << "capture         Capture records from a BUS file" << std::endl
  << "merge           Merge bus files from same experiment" << std::endl
  << "text            Convert a binary BUS file to a tab-delimited text file" << std::endl
  << "extract         Extract FASTQ reads correspnding to reads in BUS file" << std::endl
  << "linker          Remove section of barcodes in BUS files" << std::endl
  << "collapse        Turn BUS files into a BUG file" << std::endl
  << "umicorrect      Correct read errors in UMIs in BUG files" << std::endl
  << "predquant       Quantifies the counts in BUG files with unseen species prediction" << std::endl
  << "version         Prints version number" << std::endl 
  << "cite            Prints citation information" << std::endl
  << std::endl
  << "Running bustools <CMD> without arguments prints usage information for <CMD>"
  << std::endl << std::endl;
}



void Bustools_sort_Usage() {
  std::cout << "Usage: bustools sort [options] bus-files" << std::endl << std::endl
  << "Options: " << std::endl
  << "Default behavior (with no flag) is to sort by barcode, UMI, ec, then flag" << std::endl
  << "-t, --threads         Number of threads to use" << std::endl
  << "-m, --memory          Maximum memory used" << std::endl
  << "-T, --temp            Location and prefix for temporary files " << std::endl
  << "                      required if using -p, otherwise defaults to output" << std::endl 
  << "-o, --output          File for sorted output" << std::endl
  << "-p, --pipe            Write to standard output" << std::endl
  << "    --umi             Sort by UMI, barcode, then ec" << std::endl
  << "    --count           Sort by multiplicity, barcode, UMI, then ec" << std::endl
  << "    --flags           Sort by flag, barcode, UMI, then ec" << std::endl
  << std::endl;
}

void Bustools_capture_Usage() {
  std::cout << "Usage: bustools capture [options] bus-files" << std::endl << std::endl
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

void Bustools_merge_Usage() {
  std::cout << "Usage: bustools merge [options] directories" << std::endl
  << "  Note: BUS files should be sorted by flag" << std::endl << std::endl
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

void Bustools_fromtext_Usage() {
  std::cout << "Usage: bustools fromtext [options] text-files" << std::endl << std::endl
  << "Options: " << std::endl
  << "-o, --output          File for bus output" << std::endl
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
  std::cout << "Usage: bustools count [options] sorted-bus-files" << std::endl << std::endl
  << "Options: " << std::endl
  << "-o, --output          Output directory gene matrix files" << std::endl
  << "-g, --genemap         File for mapping transcripts to genes" << std::endl
  << "-e, --ecmap           File for mapping equivalence classes to transcripts" << std::endl
  << "-t, --txnames         File with names of transcripts" << std::endl
  << "    --genecounts      Aggregate counts to genes only" << std::endl
  << "    --em              Estimate gene abundances using EM algorithm" << std::endl 
  << "-m, --multimapping    Include bus records that pseudoalign to multiple genes" << std::endl
  << std::endl;
}

void Bustools_whitelist_Usage() {
  std::cout << "Usage: bustools whitelist [options] sorted-bus-file" << std::endl << std::endl
    << "Options: " << std::endl
    << "-o, --output        File for the whitelist" << std::endl
    << "-f, --threshold     Minimum number of times a barcode must appear to be included in whitelist" << std::endl
    << std::endl;
}

void Bustools_project_Usage() {
  std::cout << "Usage: bustools project [options] sorted-bus-file" << std::endl << std::endl
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

void Bustools_inspect_Usage() {
  std::cout << "Usage: bustools inspect [options] sorted-bus-file" << std::endl << std::endl
    << "Options: " << std::endl
    << "-o, --output          File for JSON output (optional)" << std::endl
    << "-e, --ecmap           File for mapping equivalence classes to transcripts" << std::endl
    << "-w, --whitelist       File of whitelisted barcodes to correct to" << std::endl
    << "-p, --pipe            Write to standard output" << std::endl
    << std::endl;
}

void Bustools_linker_Usage() {
  std::cout << "Usage: bustools linker [options] bus-files" << std::endl << std::endl
    << "Options: " << std::endl
    << "-o, --output          Output BUS file" << std::endl
    << "-s, --start           Start coordinate for section of barcode to remove (0-indexed, inclusive)" << std::endl
    << "-e, --end             End coordinate for section of barcode to remove (0-indexed, exclusive)" << std::endl
    << "-p, --pipe            Write to standard output" << std::endl
    << std::endl;
}

void Bustools_collapse_Usage() {
  std::cout << "Usage: bustools collapse [options] sorted-bus-files" << std::endl << std::endl
  << "Options: " << std::endl
  << "-o, --output          Output directory gene matrix files" << std::endl
  << "-g, --genemap         File for mapping transcripts to genes" << std::endl
  << "-e, --ecmap           File for mapping equivalence classes to transcripts" << std::endl
  << "-t, --txnames         File with names of transcripts" << std::endl
  << "-p, --pipe            Write to standard output" << std::endl
  << std::endl;
}

void Bustools_umicorrect_Usage() {
  std::cout << "Usage: bustools count [options] sorted-bus-files" << std::endl << std::endl
  << "Options: " << std::endl
  << "-o, --output          Output directory gene matrix files" << std::endl
  << "-p, --pipe            Write to standard output" << std::endl
  << std::endl;
}

void Bustools_predquant_Usage() {
  std::cout << "Usage: bustools count [options] sorted-bus-files" << std::endl << std::endl
  << "Options: " << std::endl
  << "-o, --output            Output directory gene matrix files" << std::endl
  << "-g, --genesnames        File containing gene names, output from collapse" << std::endl
  << "-m, --multimapping      Include bus records that pseudoalign to multiple genes" << std::endl
  << "-G, --goodtoulmin       Base prediction on the Good-Toulmin algorithm" << std::endl
  << "-P, --predtarget        The number of times to increase the counts for prediction" << std::endl
  << "-I, --incl_bucket_limit The number of counts needed for a gene to be included in buckets" << std::endl
  << "-U, --use_bucket_limit  Genes below this count will be predicted using buckets" << std::endl
  << "-N, --num_buckets	      Number of buckets used" << std::endl
  << std::endl;
}


void Bustools_extract_Usage() {
  std::cout << "Usage: bustools extract [options] sorted-bus-file" << std::endl
    << "  Note: BUS file should be sorted by flag using bustools sort --flag" << std::endl << std::endl
    << "Options: " << std::endl
    << "-o, --output          Output directory for FASTQ files" << std::endl
    << "-f, --fastq           FASTQ file(s) from which to extract reads (comma-separated list)" << std::endl
    << "-N, --nFastqs         Number of FASTQ file(s) per run" << std::endl
    << std::endl;
}

void print_citation() {
  std::cout << "When using this program in your research, please cite" << std::endl << std::endl
       << "  Melsted, P., Booeshaghi, A. S., et al." << std::endl
       << "  Modular and efficient pre-processing of single-cell RNA-seq, "<< std::endl
       << "  bioRxiv doi:10.1101/673285" << std::endl
       << std::endl;
}

void print_version() {
  std::cout << "bustools, version " << 	BUSTOOLS_VERSION << std::endl;
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
        //bustools_sort_orig(opt);
        bustools_sort(opt);
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
        bustools_merge(opt);
      } else {
        Bustools_merge_Usage();
        exit(1);
      }
    } else if (cmd == "cite"){
        print_citation();
    } else if (cmd == "version"){
        print_version();
    } else if (cmd == "dump" || cmd == "text") {
      if (disp_help) {
        Bustools_dump_Usage();
        exit(0);        
      }
      parse_ProgramOptions_dump(argc-1, argv+1, opt);
      if (check_ProgramOptions_dump(opt)) { //Program options are valid
	    bustools_text(opt);
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
        bustools_correct(opt);        
      } else {
        Bustools_dump_Usage();
        exit(1);
      }
    } else if (cmd == "fromtext") {
      if (disp_help) {
        Bustools_fromtext_Usage();
        exit(0);        
      }
      parse_ProgramOptions_fromtext(argc-1, argv+1, opt);
      if (check_ProgramOptions_fromtext(opt)) { //Program options are valid
	    bustools_fromtext(opt);
      } else {
        Bustools_fromtext_Usage();
        exit(1);
      }
    } else if (cmd == "count") {
      if (disp_help) {
        Bustools_count_Usage();
        exit(0);        
      }
      parse_ProgramOptions_count(argc-1, argv+1, opt);
      if (check_ProgramOptions_count(opt)) { //Program options are valid
        bustools_count(opt);
      } else {
        Bustools_count_Usage();
        exit(1);
      }
    } else if (cmd == "capture") {
      if (disp_help) {
        Bustools_capture_Usage();
        exit(0);        
      }
      parse_ProgramOptions_capture(argc-1, argv+1, opt);
      if (check_ProgramOptions_capture(opt)) { //Program options are valid
        bustools_capture(opt);
      } else {
        Bustools_capture_Usage();
        exit(1);
      }
    } else if (cmd == "whitelist") {
      if (disp_help) {
        Bustools_whitelist_Usage();
        exit(0);        
      }
      parse_ProgramOptions_whitelist(argc-1, argv+1, opt);
      if (check_ProgramOptions_whitelist(opt)) { //Program options are valid
        bustools_whitelist(opt);
      } else {
        Bustools_whitelist_Usage();
        exit(1);
      }
    } else if (cmd == "project") {
      if (disp_help) {
        Bustools_project_Usage();
        exit(0);
      }
      parse_ProgramOptions_project(argc-1, argv+1, opt);
      if (check_ProgramOptions_project(opt)) { //Program options are valid
        bustools_project(opt);
      } else {
        Bustools_project_Usage();
        exit(1);
      }
    } else if (cmd == "inspect") {
      if (disp_help) {
        Bustools_inspect_Usage();
        exit(0);
      }
      parse_ProgramOptions_inspect(argc-1, argv+1, opt);
      if (check_ProgramOptions_inspect(opt)) { //Program options are valid
        bustools_inspect(opt);
      } else {
        Bustools_inspect_Usage();
        exit(1);
      }
    } else if (cmd == "linker") {
      if (disp_help) {
        Bustools_linker_Usage();
        exit(0);
      }
      parse_ProgramOptions_linker(argc-1, argv+1, opt);
      if (check_ProgramOptions_linker(opt)) { //Program options are valid
        bustools_linker(opt);
      } else {
        Bustools_linker_Usage();
        exit(1);
      }
    } else if (cmd == "collapse") {
      if (disp_help) {
        Bustools_collapse_Usage();
        exit(0);        
      }
      parse_ProgramOptions_collapse(argc-1, argv+1, opt);
      if (check_ProgramOptions_collapse(opt)) { //Program options are valid
	    bustools_collapse(opt);
      } else {
        Bustools_collapse_Usage();
        exit(1);
      }
    } else if (cmd == "umicorrect") {
      if (disp_help) {
        Bustools_umicorrect_Usage();
        exit(0);        
      }
      parse_ProgramOptions_umicorrect(argc-1, argv+1, opt);
      if (check_ProgramOptions_umicorrect(opt)) { //Program options are valid
	    bustools_umicorrect(opt);
      } else {
        Bustools_umicorrect_Usage();
        exit(1);
      }
    } else if (cmd == "predquant") {
      if (disp_help) {
        Bustools_predquant_Usage();
        exit(0);        
      }
      parse_ProgramOptions_predquant(argc-1, argv+1, opt);
      if (check_ProgramOptions_predquant(opt)) { //Program options are valid
	    bustools_predquant(opt);
      } else {
        Bustools_predquant_Usage();
        exit(1);
      }
    } else if (cmd == "extract") {
      if (disp_help) {
        Bustools_extract_Usage();
        exit(0);
      }
      parse_ProgramOptions_extract(argc-1, argv+1, opt);
      if (check_ProgramOptions_extract(opt)) { //Program options are valid
        bustools_extract(opt);
      } else {
        Bustools_extract_Usage();
        exit(1);
      }
    } else {
      std::cerr << "Error: invalid command " << cmd << std::endl;
      Bustools_Usage();      
    }

  }
}
