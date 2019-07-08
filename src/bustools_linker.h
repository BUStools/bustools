#include "Common.hpp"

typedef struct lnk_Record {
  uint64_t barcode;
  uint64_t UMI;
  uint32_t ec;
} lnk_Record;

void bustools_linker(Bustools_opt &opt);
