#include "Common.hpp"

typedef struct wl_Record {
  uint64_t barcode;
  uint32_t R;
  uint32_t U;
  uint32_t count;

  wl_Record(uint64_t b, uint32_t r, uint32_t u, uint32_t c) :
    barcode(b), R(r), U(u), count(c) {}
} wl_Record;

void bustools_whitelist(Bustools_opt &opt);
