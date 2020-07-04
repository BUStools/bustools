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

bool parse_ProjectMap(const std::string &filename, std::unordered_map<uint64_t, uint64_t> &project_map);

uint64_t stringToBinary(const std::string &s, uint32_t &flag);
uint64_t stringToBinary(const char* s, const size_t len, uint32_t &flag);
std::string binaryToString(uint64_t x, size_t len);
int hamming(uint64_t a, uint64_t b, size_t len);

void write_ecs_block(std::ostream& o, const std::vector<int32_t>& v);
void write_bug_entry(std::ostream& o, const BUSData& bd, const std::vector<int32_t>* pGeneList = nullptr);
void add_ecs_from_block(const BUSData& bd, std::vector<int32_t>& v, size_t& leftToRead);

bool parseBugGenes(const std::string& filename, std::vector<std::string>& genes);
/*
class BUGReader {
public:
	BUGReader(const std::vector<std::string>& filenames, bool bUseStdIn) : m_filenames(filenames), m_bUseStdIn(bUseStdIn){
		openNextFile();
	}

	//For performance reasons (to avoid copying data), this class returns a pointer to the data
	//it only uses the vector in case of multiple genes, so the vector should only be used if the
	//ec field is negative
	BUSData* getBug(std::vector<int32_t>& genes) {

	}

private:
	size_t m_nextFile = 0;
	size_t pos = 0;
	std::vector<std::string> m_filenames;
	bool m_bUseStdIn = false;
	std::shared_ptr<std::istream> m_pIn;
	std::shared_ptr<std::ifstream> m_pInf;

	bool openNextFile() {
		std::streambuf* inbuf;
		if (!m_bUseStdIn) {
			if (m_nextFile >= m_filenames.size()) {
				return false;
			}
			m_pInf.reset(new std::ifstream());
			m_pInf->open(m_filenames[m_nextFile].c_str(), std::ios::binary);
			inbuf = m_pInf->rdbuf();
			if (!*m_pInf) {
				std::cerr << "Failed to read file: " << m_filenames[m_nextFile].c_str() << std::endl;
			}
		}
		else {
			if (m_nextFile > 0) {
				return false;
			}
			inbuf = std::cin.rdbuf();
		}
		m_pIn.reset(new std::istream(inbuf));
		m_nextFile++;
		return true;
	}

};
*/
#endif // KALLISTO_BUSDATA_H
