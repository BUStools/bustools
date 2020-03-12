#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <queue>
#include <functional>

#include "Common.hpp"
#include "BUSData.h"
#include "bustools_text.h"


void write_bus(std::ostream& o, std::string& bc, const std::string& umi, const std::string& countsstr, const std::string& ecstr) {
	static uint32_t f = 0;
	static std::vector<int32_t> ecs;
	BUSData bd;
	bd.barcode = stringToBinary(bc, f);
	bd.UMI = stringToBinary(umi, f);
	bd.count = atoi(countsstr.c_str());

	//check for commas in ecstr
	split_string(ecstr, ecs);
	write_bug_entry(o, bd, &ecs);
}

void bustools_text(const Bustools_opt& opt) {
	if (!opt.text_import) {
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
	} else {
		BUSHeader h;

		std::streambuf* buf = nullptr;
		std::ofstream of;
		int nr = 0;

		if (!opt.stream_out) {
			of.open(opt.output, std::ios::binary);
			buf = of.rdbuf();
		}
		else {
			buf = std::cout.rdbuf();
		}
		std::ostream o(buf);

		std::string bc, umi, countsstr, ecstr;

		bool bHeaderWritten = false;


		//	char magic[4];
		uint32_t version = 0;
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

			if (!bHeaderWritten) {
				bHeaderWritten = true;
				//This is a trick to figure out the length of the barcode and umi in the header
				in >> bc >> umi >> ecstr >> countsstr;
				h.bclen = uint32_t(bc.length());
				h.umilen = uint32_t(umi.length());
				h.version = BUSFORMAT_VERSION;
				writeHeader(o, h);
				//don't forget to write the item just read
				write_bus(o, bc, umi, countsstr, ecstr);
				++nr;
			}
			int rc = 0;
			while (in >> bc >> umi >> ecstr >> countsstr) {
				write_bus(o, bc, umi, countsstr, ecstr);
				++nr;
			}
		}
		if (!opt.stream_out) {
			of.close();
		}
		std::cerr << "Read in " << nr << " text records" << std::endl;

	}
}
