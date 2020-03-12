#ifndef BUSTOOLS_COMMON_HPP
#define BUSTOOLS_COMMON_HPP

#include <cassert>
#include <cmath>
#include <algorithm>
#include <stdint.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <sstream>


#define BUSTOOLS_VERSION "0.39.4"

enum CAPTURE_TYPE : char {CAPTURE_NONE = 0, CAPTURE_TX, CAPTURE_BC, CAPTURE_UMI, CAPTURE_F};
enum SORT_TYPE : char {SORT_BC = 0, SORT_UMI, SORT_F, SORT_COUNT};
enum PROJECT_TYPE : char { PROJECT_BC = 0, PROJECT_UMI, PROJECT_TX, PROJECT_F };
enum class PREDQUANT_ALG { GOOD_TOULMIN, PRESEQ };

struct Bustools_opt {
  int threads;
  
  std::string whitelist;  
  std::string output;
  std::vector<std::string> files;

  bool stream_in = false;
  bool stream_out = false;

  /* extract */
  int nFastqs;
  std::vector<std::string> fastq;
  
  char type;

  int ec_d;
  int ec_dmin;
  size_t max_memory;
  std::string temp_files;

  /* count, and other things */
  std::string count_genes;
  std::string count_ecs;
  std::string count_txp;
  bool count_em = false;
  bool count_collapse = false;
  bool count_gene_multimapping = false;

  /* project */
  std::string map;
  std::string output_folder;
  
  /* capture */
  std::string capture;
  bool complement = false;
  bool filter = false;

  /* whitelist */
  int threshold;

  /* text */
  bool text_dumpflags = false;

  /* linker */
  int start, end;

  // predquant
  std::string predquant_genefile;
  PREDQUANT_ALG predquant_algorithm = PREDQUANT_ALG::PRESEQ;
  int predquant_pred_target = 10;
  int predquant_num_buckets = 100;
  double predquant_use_bucket_limit = 200; //if the prediction should be mapped to a bucket or not (number of UMIs in total per gene)
  double predquant_include_bucket_limit = 200;//if gene should be used in buckets for other genes or not (number of UMIs in total per gene)

  Bustools_opt() : threads(1), max_memory(1ULL<<32), type(0),
    threshold(0), start(-1), end(-1)  {}
};

static const char alpha[4] = {'A','C','G','T'};

inline size_t rndup(size_t v) {

  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v |= v >> 32;
  v++;

  return v;
}

inline uint32_t rndup(uint32_t v) {

  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;

  return v;
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
std::vector<int32_t> intersect(std::vector<int32_t> &u, std::vector<int32_t> &v);
std::vector<int32_t> union_vectors(const std::vector<std::vector<int32_t>> &v);
std::vector<int32_t> intersect_vectors(const std::vector<std::vector<int32_t>> &v);
int32_t intersect_ecs(const std::vector<int32_t> &ecs, std::vector<int32_t> &u, const std::vector<int32_t> &genemap, std::vector<std::vector<int32_t>> &ecmap, std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> &ecmapinv, std::vector<std::vector<int32_t>> &ec2genes);
void vt2gene(const std::vector<int32_t> &v, const std::vector<int32_t> &genemap, std::vector<int32_t> &glist);
void intersect_genes_of_ecs(const std::vector<int32_t> &ecs, const  std::vector<std::vector<int32_t>> &ec2genes, std::vector<int32_t> &glist);
int32_t intersect_ecs_with_genes(const std::vector<int32_t> &ecs, const std::vector<int32_t> &genemap, std::vector<std::vector<int32_t>> &ecmap, std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> &ecmapinv, std::vector<std::vector<int32_t>> &ec2genes, bool assumeIntersectionIsEmpty = true);
void create_ec2genes(const std::vector<std::vector<int32_t>> &ecmap, const std::vector<int32_t> &genemap, std::vector<std::vector<int32_t>> &ec2gene);

void split_string(const std::string& str, std::vector<int32_t>& res);


#endif // BUSTOOLS_COMMON_HPP
