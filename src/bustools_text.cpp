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


	//	char magic[4];
	uint32_t version = 0;
	for (const auto& infn : opt.files) {
		std::streambuf* inbuf;
		std::ifstream inf;
		if (!opt.stream_in) {
			inf.open(infn.c_str(), std::ios::binary);
			inbuf = inf.rdbuf();
			if (!inf) {
				std::cerr << "Failed to read file: " << infn << std::endl;
			}
		}
		else {
			inbuf = std::cin.rdbuf();
		}
		std::istream in(inbuf);


		parseHeader(in, h);
		uint32_t bclen = h.bclen;
		uint32_t umilen = h.umilen;
		std::vector<int32_t> ecs;
		BUSData bdmult;
		int rc = 0;
		size_t leftToRead = 0;
		while (true) {
			in.read((char*)p, N * sizeof(BUSData));
			size_t rc = in.gcount() / sizeof(BUSData);
			if (rc == 0) {
				break;
			}
			nr += rc;
			for (size_t i = 0; i < rc; i++) {
				if (leftToRead > 0) {
					add_ecs_from_block(p[i], ecs, leftToRead);
					if (leftToRead == 0) {
						//write as text
						o << binaryToString(bdmult.barcode, bclen) << "\t" << binaryToString(bdmult.UMI, umilen) << "\t";
						bool bSkipComma = true;
						for (size_t e = 0; e < ecs.size(); ++e) {
							if (!bSkipComma) {
								o << ',';
							}
							bSkipComma = false;
							o << ecs[e];
						}
							
						o << "\t" << bdmult.count << "\n";
					}
				} else {
					if (p[i].ec < 0) {
						//start a "multiblock" read
						ecs.clear();
						bdmult = p[i];
						leftToRead = -p[i].ec;
					} else {
						o << binaryToString(p[i].barcode, bclen) << "\t" << binaryToString(p[i].UMI, umilen) << "\t" << p[i].ec << "\t" << p[i].count;
						if (opt.text_dumpflags) {
							o << "\t" << p[i].flags;
						}
						o << "\n";
					}
				}
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
	std::string line, bc, umi, ecstr;
	int32_t count;
	std::vector<int32_t> ecs;

	for (const auto& infn : opt.files) {
		std::streambuf* inbuf;
		std::ifstream inf;
		if (!opt.stream_in) {
			inf.open(infn.c_str()); //do not open as binary, not sure what effect it will have in practice though
			inbuf = inf.rdbuf();
		}
		else {
			inbuf = std::cin.rdbuf();
		}
		std::istream in(inbuf);
		
		while (std::getline(in, line)) {
			std::stringstream ss(line);
			ss >> bc >> umi >> ecstr >> count;
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
			b.count = count;
			b.flags = 0;
			//check for commas in ecstr
			split_string(ecstr, ecs);
			write_bug_entry(o, b, &ecs);
			++nr;
		}
	}
	if (!opt.stream_out) {
		of.close();
	}
	std::cerr << "Read in " << nr << " text records" << std::endl;
}
