#include "BUSData.h"


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


bool parseHeader(std::ifstream &inf, BUSHeader &header) {
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

bool writeHeader(std::ofstream &outf, const BUSHeader &header) {
  outf.write("BUS\0", 4);
  outf.write((char*)(&header.version), sizeof(header.version));
  outf.write((char*)(&header.bclen), sizeof(header.bclen));
  outf.write((char*)(&header.umilen), sizeof(header.umilen));
  uint32_t tlen = header.text.size();
  outf.write((char*)(&tlen), sizeof(tlen));
  outf.write((char*)header.text.c_str(), tlen);

  return true;
}