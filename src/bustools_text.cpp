#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <queue>
#include <functional>

#include "Common.hpp"
#include "BUSData.h"
#include "bustools_text.h"


void bustools_text(const Bustools_opt& opt) {
	BUSHeader h;
	size_t nr = 0;
	size_t N = 100000;
	BUSData* p = new BUSData[N];

	std::streambuf* buf = nullptr;
	std::ofstream of;

	if (!opt.stream_out) {
		of.open(opt.output);
		buf = of.rdbuf();
	}
	else {
		buf = std::cout.rdbuf();
	}
	std::ostream o(buf);


	uint32_t version = 0;
	for (const auto& infn : opt.files) {
		std::streambuf* inbuf;
		std::ifstream inf;
		if (!opt.stream_in) {
			inf.open(infn.c_str(), std::ios::binary);
			inbuf = inf.rdbuf();
		}
		else {
			inbuf = std::cin.rdbuf();
		}
		std::istream in(inbuf);


		parseHeader(in, h);
		uint32_t bclen = h.bclen;
		uint32_t umilen = h.umilen;
		int rc = 0;
		while (true) {
			in.read((char*)p, N * sizeof(BUSData));
			size_t rc = in.gcount() / sizeof(BUSData);
			if (rc == 0) {
				break;
			}
			nr += rc;
			for (size_t i = 0; i < rc; i++) {
				o << binaryToString(p[i].barcode, bclen) << "\t" << binaryToString(p[i].UMI, umilen) << "\t" << p[i].ec << "\t" << p[i].count;
				if (opt.text_dumpflags) {
					o << "\t" << p[i].flags;
				}
				if (opt.text_dumppad)
				{
				  o << "\t" << p[i].pad;
				}
				o << "\n";
			}
		}
	}
	delete[] p; p = nullptr;
	if (!opt.stream_out) {
		of.close();
	}
	std::cerr << "Read in " << nr << " BUS records" << std::endl;
}

void bustools_fromtext(const Bustools_opt& opt) {
	std::streambuf* buf = nullptr;
	std::ofstream of;
	int nr = 0;

	if (!opt.stream_out) {
		of.open(opt.output, std::ios::binary);
		buf = of.rdbuf();
	} else {
		buf = std::cout.rdbuf();
	}
	std::ostream o(buf);

	BUSHeader h;
	uint32_t f;
	bool out_header_written = false;
	std::string line, bc, umi;
	int32_t ec = 0, count = 0;

	for (const auto& infn : opt.files) {
		std::streambuf* inbuf;
		std::ifstream inf;
		if (!opt.stream_in) {
			inf.open(infn.c_str()); //do not open as binary, will give problems with getline in Windows
			inbuf = inf.rdbuf();
		}
		else {
			inbuf = std::cin.rdbuf();
		}
		std::istream in(inbuf);
		
		while (std::getline(in, line)) {
			//allow for comments and empty rows:
			if (line.empty() || line[0] == '#') {
				continue;
			}
			std::stringstream ss(line);
			//this will automatically allow for comments after count on each line, as long as there is a whitespace in between
			if (ss >> bc >> umi >> ec >> count) { //if the if is not here, empty lines at the end  of the file will add extra entries identical to the last one...
				if (!out_header_written) {
					h.bclen = (uint32_t)bc.size();
					h.umilen = (uint32_t)umi.size();
					h.version = BUSFORMAT_VERSION;
					h.text = "converted from text format";
					writeHeader(o, h);
					out_header_written = true;
				}
				BUSData b;
				b.barcode = stringToBinary(bc, f);
				b.UMI = stringToBinary(umi, f);
				b.ec = ec;
				b.count = count;
				b.flags = 0;
				o.write((char*)&b, sizeof(b));
				++nr;
			}
		}
	}
	if (!opt.stream_out) {
		of.close();
	}
	std::cerr << "Read in " << nr << " text records" << std::endl;
}
