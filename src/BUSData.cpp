#include "BUSData.h"

#include <cstring>
#include <assert.h>
#include <unordered_map>
#include <sstream>
#include <iostream>
#include <algorithm>

uint64_t stringToBinary(const std::string &s, uint32_t &flag) {
  return stringToBinary(s.c_str(), s.size(), flag);
}

std::string binaryToString(uint64_t x, size_t len) {
  std::string s(len, 'N');
  size_t sh = len-1;
  for (size_t i = 0; i < len; i++) {
    char c = 'N';
    switch((x >> (2*sh)) & 0x03ULL) {
      case 0x00: c = 'A'; break;
      case 0x01: c = 'C'; break;
      case 0x02: c = 'G'; break;
      case 0x03: c = 'T'; break;
    }
    sh--;
    s.at(i) = c;
  }
  return std::move(s);
}


int hamming(uint64_t a, uint64_t b, size_t len) {
  uint64_t df = a^b;
  int d = 0;
  size_t sh = len-1;
  for (size_t i = 0; i < len; ++i) {
    if (((df>>(2*sh)) & 0x03ULL) != 0) {
      d++;
    }    
    sh--;
  }
  return d;
}

uint64_t stringToBinary(const char* s, const size_t len, uint32_t &flag) {
  uint64_t r = 0;
  flag = 0;
  int numN = 0;
  size_t posN = 0;
  size_t k = len;
  if (k > 32) {
    k = 32;
  }
  for (size_t i = 0; i < k; ++i) {
    uint64_t x = ((*s) & 4) >> 1;
    if (((*s) & 3) == 2) {
      if (numN == 0) {
        posN = i;
      }
      ++numN;
    }
    r = r << 2;
    r |= (x + ((x ^ (*s & 2)) >>1));
    s++;
  }
  if (numN>0) {
    if (numN > 3) {
      numN = 3;      
    }    
    flag = (numN & 3) | (posN & 15) << 2;
  }
  return r;
}


bool parseHeader(std::istream &inf, BUSHeader &header) {
  char magic[4];  
  inf.read((char*)(&magic[0]), 4);
  if (std::strcmp(&magic[0], "BUS\0") != 0) {
    return false;
  }
  inf.read((char*)(&header.version), sizeof(header.version));
  if (header.version != BUSFORMAT_VERSION) {
    return false;
  }
  inf.read((char*)(&header.bclen), sizeof(header.bclen));
  inf.read((char*)(&header.umilen), sizeof(header.umilen));
  uint32_t tlen = 0;
  inf.read((char*)(&tlen), sizeof(tlen));
  char* t = new char[tlen+1];
  inf.read(t, tlen);
  t[tlen] = '\0';
  header.text.assign(t);
  delete[] t;

  return true;
}



bool parseECs(const std::string &filename, BUSHeader &header) {
  auto &ecs = header.ecs; 
  std::ifstream inf(filename.c_str());
  std::string line, t;
  line.reserve(10000);
  
  std::vector<int32_t> c;
  
  int i = 0;
  while (std::getline(inf, line)) {       
    c.clear();
    int ec = -1;
    if (line.size() == 0) {
      continue;
    }
    std::stringstream ss(line);
    ss >> ec;
    assert(ec == i);
    while (std::getline(ss, t, ',')) {
      c.push_back(std::stoi(t));
    }

    ecs.push_back(std::move(c));
    i++;
  }

  return true;
}

bool writeECs(const std::string &filename, const BUSHeader &header) {
  std::ofstream outf;
  outf.open(filename.c_str(), std::ios::out);

  if (!outf.is_open()) {
    return false;
  }

  
  size_t n = header.ecs.size();
  for (size_t ec = 0; ec < n; ec++) {
    const auto& v = header.ecs[ec];
    outf << ec << "\t";
    bool first = true;
    for (auto x : v) {
      if (!first) {
        outf << ",";
      } else {
        first = false;
      }
      outf << x;
    }
    outf << "\n";
  }
  outf.close();

  return true;
}

bool writeGenes(const std::string &filename, const std::unordered_map<std::string, int32_t>  &genenames) {
  std::ofstream outf;
  outf.open(filename.c_str(), std::ios::out);

  if (!outf.is_open()) {
    return false;
  }
  std::vector<std::string> names;
  names.resize(genenames.size());
  for (const auto &x : genenames) {
    if (x.second >= 0) {
      names[x.second] = x.first;
    }
  }
  for (const auto &x : names) {
    outf << x << "\n";
  }

  return true;
}

bool parseTranscripts(const std::string &filename, std::unordered_map<std::string, int32_t> &txnames) {
  std::ifstream inf(filename.c_str());

  int i = 0;
  std::string txp;
  while (inf >> txp) {
    txnames.insert({txp, i});
    i++;
  }
  return true;
}

bool parseTxCaptureList(const std::string &filename, std::unordered_map<std::string, int32_t> &txnames, std::unordered_set<uint64_t> &captures) {
  std::ifstream inf(filename.c_str());

  std::string txp;
  while (inf >> txp) {
    auto it = txnames.find(txp);
    if (it == txnames.end()) {
      std::cerr << "Error: could not find capture transcript " << txp << " in transcript list" << std::endl;
      return false;
    } 
    captures.insert((uint64_t) it->second);
  }
  return true;
}

bool parseBcUmiCaptureList(const std::string &filename, std::unordered_set<uint64_t> &captures) {
  std::ifstream inf(filename.c_str());

  std::string inp;
  uint32_t flag; // Unused
  while (getline(inf, inp)) {
    captures.insert(stringToBinary(inp, flag));
  }

  return true;
}

bool parse_ProjectMap(const std::string &filename, std::unordered_map<uint64_t, uint64_t> &project_map) {
  // This function occurs in 3 places: here, BUSData.h, and bustools_project.cpp
  std::ifstream inf(filename.c_str());

  std::string line, t; // whats the point of the t?
  line.reserve(10000);
  uint32_t flag; // unused
  while (std::getline(inf, line)) {
    std::stringstream ss(line);
    std::string source, dest;

    ss >> source >> dest; 
    project_map[stringToBinary(source, flag)] = stringToBinary(dest, flag);
  }
  return true;
}

bool parseFlagsCaptureList(const std::string &filename, std::unordered_set<uint64_t> &captures) {
  std::ifstream inf(filename.c_str());
  
  std::string inp;
  while (getline(inf, inp)) {
    captures.insert(stoi(inp));
  }

  return true;
}

bool parseGenes(const std::string &filename, const std::unordered_map<std::string, int32_t> &txnames, std::vector<int32_t> &genemap, std::unordered_map<std::string, int32_t> &genenames) {
  std::ifstream inf(filename.c_str());

  std::string line, t;
  line.reserve(10000);

  int i = 0;
  while (std::getline(inf,line)) {
    std::stringstream ss(line);
    std::string txp, gene;
    ss >> txp >> gene;
    auto it = txnames.find(txp);
    if (it != txnames.end()) {
      auto i = it->second;
      auto git = genenames.find(gene);
      auto gi = -1;
      if (git == genenames.end()) {
        gi = genenames.size();
        genenames.insert({gene,gi});
      } else {
        gi = git->second;
      }
      genemap[i] = gi;
    }
  }

  return true;
}

bool writeHeader(std::ostream &outf, const BUSHeader &header) {
  outf.write("BUS\0", 4);
  outf.write((char*)(&header.version), sizeof(header.version));
  outf.write((char*)(&header.bclen), sizeof(header.bclen));
  outf.write((char*)(&header.umilen), sizeof(header.umilen));
  uint32_t tlen = header.text.size();
  outf.write((char*)(&tlen), sizeof(tlen));
  outf.write((char*)header.text.c_str(), tlen);

  return true;
}

void write_ecs_block(std::ostream& o, const std::vector<int32_t>& v)
{
	static char zeros[32] = { 0 };
	o.write((char*)&v[0], sizeof(int32_t) * v.size());
	const size_t numsPerBug = sizeof(BUSData) / sizeof(int32_t);
	size_t numodd = v.size() % numsPerBug;
	if (numodd > 0) {
		size_t padding = numsPerBug - numodd;
		o.write(zeros, sizeof(int32_t) * padding);
	}
}

//If pGeneList is not null, the value in bd.ec will be ignored and controlled by the gene list.
//Otherwise, the bd.ec value will be written as is
void write_bug_entry(std::ostream& o, const BUSData& bd, const std::vector<int32_t>* pGeneList)
{
	if (pGeneList) {
		auto bdc = bd;//must copy it to be able to modify it. Direct modification would be unexpected.
		auto ecs = *pGeneList;
		//if no commas, just write ec, otherwise write negative number of ecs and thereafter a list
		if (ecs.size() == 1) {
			bdc.ec = *ecs.begin();
			o.write((char*)&bdc, sizeof(bd));
		} else {
			//This is where it gets a little more complicated. If we have multiple genes,
			//we write the negative number of genes (i.e. if 3 genes, we write -3)
			//We then write blocks of the same size as the BUS entry, as many as we need

			bdc.ec = -int32_t(ecs.size());
			o.write((char*)&bdc, sizeof(bd));

			write_ecs_block(o, ecs);
		}
	} else {
		o.write((char*)&bd, sizeof(bd));
	}
}


void add_ecs_from_block(const BUSData& bd, std::vector<int32_t>& v, size_t& leftToRead) {
	const size_t numsPerBug = sizeof(BUSData) / sizeof(int32_t);
	auto numToGrab = std::min(leftToRead, numsPerBug);
	leftToRead -= numToGrab;
	auto p = (int32_t*)(&bd);
	for (size_t i = 0; i < numToGrab; ++i) {
		v.push_back(p[i]);
	}
}

bool parseBugGenes(const std::string& filename, std::vector<std::string>& genes)
{
	std::ifstream inf(filename.c_str());
	if (!inf) {
		std::cerr << "Failed to load gene file: " << filename << std::endl;
		return false;
	}

	std::string line, t;
	line.reserve(10000);

	int i = 0;
	while (std::getline(inf, line)) {
		genes.push_back(line);
	}

	return true;
}
