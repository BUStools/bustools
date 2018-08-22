#ifndef BUSTOOLS_COMMON_HPP
#define BUSTOOLS_COMMON_HPP

#include <cassert>
#include <algorithm>
#include <stdint.h>
#include <vector>
#include <string>


#define BUSTOOLS_VERSION "0.1"


struct Bustools_opt {
  int threads;

  std::vector<std::string> files;

  Bustools_opt() : threads(1) {}
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



#endif // BUSTOOLS_COMMON_HPP
