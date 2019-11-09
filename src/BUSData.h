#ifndef KALLISTO_BUSDATA_H
#define KALLISTO_BUSDATA_H

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stdint.h>
#include <fstream>

const uint32_t BUSFORMAT_VERSION = 1;

struct BUSTranscript {
  std::string name;
  uint32_t transcriptLength;
  BUSTranscript() : transcriptLength(0) {}
  BUSTranscript(std::string n) : name(n), transcriptLength(0) {}
};


struct BUSHeader {
  std::string text;
  std::vector<BUSTranscript> transcripts;
  std::vector<std::vector<int32_t>> ecs;
  uint32_t version;
  uint32_t bclen;
  uint32_t umilen;
  BUSHeader() : version(0), bclen(0), umilen(0) {}
};

struct BUSData {
  uint64_t barcode;
  uint64_t UMI;
  int32_t ec;
  uint32_t count;
  uint32_t flags;
  uint32_t pad;
  BUSData() : barcode(0), UMI(0), ec(-1), count(0), flags(0), pad(0) {}
};


bool parseHeader(std::istream &inf, BUSHeader &header);
bool writeHeader(std::ostream &outf, const BUSHeader &header);


bool parseECs(const std::string &filename, BUSHeader &header);
bool writeECs(const std::string &filename, const BUSHeader &header);
bool writeGenes(const std::string &filename, const std::unordered_map<std::string, int32_t>  &genenames);
bool parseGenes(const std::string &filename, const std::unordered_map<std::string, int32_t> &txnames, std::vector<int32_t> &genemap, std::unordered_map<std::string, int32_t> &genenames);
bool parseTxCaptureList(const std::string &filename, std::unordered_map<std::string, int32_t> &txnames, std::unordered_set<uint64_t> &captures);
bool parseBcUmiCaptureList(const std::string &filename, std::unordered_set<uint64_t> &captures);
bool parseFlagsCaptureList(const std::string &filename, std::unordered_set<uint64_t> &captures);
bool parseTranscripts(const std::string &filename, std::unordered_map<std::string, int32_t> &txnames);

uint64_t stringToBinary(const std::string &s, uint32_t &flag);
uint64_t stringToBinary(const char* s, const size_t len, uint32_t &flag);
std::string binaryToString(uint64_t x, size_t len);
int hamming(uint64_t a, uint64_t b, size_t len);
#endif // KALLISTO_BUSDATA_H
