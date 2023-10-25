#include "BUSData.h"

#include <cstring>
#include <assert.h>
#include <unordered_map>
#include <sstream>
#include <iostream>

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

int identifyParseHeader(std::istream &inf, BUSHeader &header, compressed_BUSHeader &comp_header)
{
  int ret = -1;
  char magic[5];

  magic[4] = '\0';

  inf.read(magic, 4);
  for (int i = 0; i < 4; ++i)
  {
    inf.putback(magic[3 - i]);
  }

  if (std::strcmp(magic, "BUS\0") == 0)
  {
    return BUSFILE_TYPE::BUSFILE * parseHeader(inf, header);
  }
  else if (std::strcmp(magic, "BUS\1") == 0)
  {
    return BUSFILE_TYPE::BUSFILE_COMPRESED * parseCompressedHeader(inf, comp_header);
  }
  else if(std::strcmp(magic, "BZI\0") == 0){
    return BUSFILE_TYPE::BUSZ_INDEX;
  }
  else if (std::strcmp(magic, "BEC\0") == 0)
  {
    return BUSFILE_TYPE::EC_MATRIX_COMPRESSED;
  }
  return BUSFILE_TYPE::EC_MATRIX;
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

bool parseCompressedHeader(std::istream &inf, compressed_BUSHeader &compheader)
{
  char magic[5];
  magic[4] = '\0';

  BUSHeader &header = compheader.bus_header;
  inf.read(magic, 4);
  if (std::strcmp(magic, "BUS\1") != 0)
  {
    std::cerr << "Invalid header magic\n";
    return false;
  }
  inf.read((char *)(&header.version), sizeof(header.version));
  if (header.version != BUSFORMAT_VERSION)
  {
    return false;
  }
  inf.read((char *)(&header.bclen), sizeof(header.bclen));
  inf.read((char *)(&header.umilen), sizeof(header.umilen));
  uint32_t tlen = 0;
  inf.read((char *)(&tlen), sizeof(tlen));
  char *t = new char[tlen + 1];
  inf.read(t, tlen);
  t[tlen] = '\0';
  header.text.assign(t);
  delete[] t;

  // We store the compressed_header-specific information after the regular header
  inf.read((char *)&compheader.chunk_size, sizeof(compheader.chunk_size));
  inf.read((char *)&compheader.pfd_blocksize, sizeof(compheader.pfd_blocksize));
  inf.read((char *)&compheader.lossy_umi, sizeof(compheader.lossy_umi));

  return true;
}

bool parseECs_stream(std::istream &in, BUSHeader &header)
{
  auto &ecs = header.ecs;
  std::string line, t;
  line.reserve(10000);

  std::vector<int32_t> c;

  int i = 0;
  bool has_reached = false;
  while (std::getline(in, line))
  {
    c.clear();
    int ec = -1;
    if (line.size() == 0)
    {
      continue;
    }
    std::stringstream ss(line);
    ss >> ec;
    assert(ec == i);
    while (std::getline(ss, t, ','))
    {
      c.push_back(std::stoi(t));
    }
    if (!has_reached)
    {
      has_reached |= !(c.size() == 1 && c[0] == i);
      if (has_reached)
      {
        std::cerr << "first line is " << i << '\n';
      }
    }

    ecs.push_back(std::move(c));
    ++i;
  }
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

bool writeGenes(const std::string &filename, const u_map_<std::string, int32_t>  &genenames) {
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

bool parseTranscripts(const std::string &filename, u_map_<std::string, int32_t> &txnames) {
  std::ifstream inf(filename.c_str());

  int i = 0;
  std::string txp;
  while (inf >> txp) {
    txnames.insert({txp, i});
    i++;
  }
  return true;
}

bool parseTxCaptureList(const std::string &filename, u_map_<std::string, int32_t> &txnames, std::unordered_set<uint64_t> &captures) {
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
    uint64_t binary_val;
    if (inp.back() == '*') {
	    inf.pop_back(); // Remove * from end of string
	    binary_val = stringToBinary(inp, flag);
	    binary_val |= (static_cast<uint64_t>(1) << 63); // Set MSB to 1 if * found at end of string (to label the barcode as "prefix")
    } else {
	    binary_val = stringToBinary(inp, flag);
    }
    captures.insert(binary_val);
  }

  return true;
}

bool parse_ProjectMap(const std::string &filename, u_map_<uint64_t, uint64_t> &project_map) {
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

bool parseGenes(const std::string &filename, const u_map_<std::string, int32_t> &txnames, std::vector<int32_t> &genemap, u_map_<std::string, int32_t> &genenames) {
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

bool parseGenesList(const std::string& filename, std::vector<std::string>& geneNames) {
	std::ifstream inf(filename.c_str());
	geneNames.clear();
	geneNames.reserve(10000);
	std::string str;
	while (inf >> str) {
		geneNames.push_back(str);
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

bool writeCompressedHeader(std::ostream &outf, const compressed_BUSHeader &compheader)
{
  outf.write("BUS\1", 4);

  // We start writing out the contents of the general header
  const auto header = compheader.bus_header;
  outf.write((char *)(&header.version), sizeof(header.version));
  outf.write((char *)(&header.bclen), sizeof(header.bclen));
  outf.write((char *)(&header.umilen), sizeof(header.umilen));
  uint32_t tlen = header.text.size();
  outf.write((char *)(&tlen), sizeof(tlen));
  outf.write((char *)header.text.c_str(), tlen);

  // We end by writing out the compressed-header-specific data
  outf.write((char *)(&compheader.chunk_size), sizeof(compheader.chunk_size));
  outf.write((char *)&compheader.pfd_blocksize, sizeof(compheader.pfd_blocksize));
  outf.write((char *)(&compheader.lossy_umi), sizeof(compheader.lossy_umi));

  return true;
}
