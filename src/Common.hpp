#ifndef BUSTOOLS_COMMON_HPP
#define BUSTOOLS_COMMON_HPP

#include <cassert>
#include <algorithm>
#include <stdint.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <sstream>


#define BUSTOOLS_VERSION "0.39.2"


struct Bustools_opt {
  int threads;
  std::string ecf;
  std::string output;
  std::string whitelist;  
  std::vector<std::string> files;

  int ec_d;
  int ec_dmin;
  size_t max_memory;
  std::string temp_files;

  std::string count_genes;
  std::string count_ecs;
  std::string count_txp;
  bool count_collapse = false;
  bool count_gene_multimapping = false;

  std::string capture;

  bool stream_in = false;
  bool stream_out = false;

  int threshold;

  Bustools_opt() : threads(1), max_memory(1ULL<<32) {}
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
int32_t intersect_ecs(const std::vector<int32_t> &ecs, std::vector<int32_t> &u, std::vector<std::vector<int32_t>> &ecmap, std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> &ecmapinv);
void vt2gene(const std::vector<int32_t> &v, const std::vector<int32_t> &genemap, std::vector<int32_t> &glist);
void intersect_genes_of_ecs(const std::vector<int32_t> &ecs, const  std::vector<std::vector<int32_t>> &ec2genes, std::vector<int32_t> &glist);
int32_t intersect_ecs_with_genes(const std::vector<int32_t> &ecs, const std::vector<int32_t> &genemap, std::vector<std::vector<int32_t>> &ecmap, std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> &ecmapinv, std::vector<std::vector<int32_t>> &ec2genes, bool assumeIntersectionIsEmpty = true);
void create_ec2genes(const std::vector<std::vector<int32_t>> &ecmap, const std::vector<int32_t> &genemap, std::vector<std::vector<int32_t>> &ec2gene);




#endif // BUSTOOLS_COMMON_HPP
