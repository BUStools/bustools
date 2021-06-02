#include <iostream>
#include <fstream>
#include <cstring>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_count.h"
#include <random>


void bustools_count(Bustools_opt &opt) {
  BUSHeader h;
  size_t nr = 0;
  size_t N = 100000;
  uint32_t bclen = 0;
  BUSData* p = new BUSData[N];

  // read and parse the equivalence class files

  std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> ecmapinv;
  std::vector<std::vector<int32_t>> ecmap;

  std::unordered_map<std::string, int32_t> txnames;
  parseTranscripts(opt.count_txp, txnames);
  std::vector<int32_t> genemap(txnames.size(), -1);
  std::unordered_map<std::string, int32_t> genenames;
  parseGenes(opt.count_genes, txnames, genemap, genenames);
  parseECs(opt.count_ecs, h);
  ecmap = std::move(h.ecs);
  ecmapinv.reserve(ecmap.size());
  for (int32_t ec = 0; ec < ecmap.size(); ec++) {
    ecmapinv.insert({ecmap[ec], ec});
  }
  std::vector<std::vector<int32_t>> ec2genes;        
  create_ec2genes(ecmap, genemap, ec2genes);


  std::ofstream of;
  std::string mtx_ofn = opt.output + ".mtx";
  std::string barcodes_ofn = opt.output + ".barcodes.txt";
  std::string ec_ofn = opt.output + ".ec.txt";
  std::string gene_ofn = opt.output + ".genes.txt";
  std::string hist_ofn = opt.output + ".hist.txt";
  std::string cu_per_cell_ofn = opt.output + ".CUPerCell.txt";
  std::string cu_ofn = opt.output + ".cu.txt";
  of.open(mtx_ofn);

  // write out the initial header
  // keep the number of newlines constant, this way it will work for both Windows and Linux
  std::string headerBuf;
  std::stringstream ssHeader(headerBuf);
  std::string headerComments = "%%MatrixMarket matrix coordinate real general\n%\n";
  ssHeader << headerComments;
  ssHeader << std::string(66, '%') << '\n';
  size_t headerLength = ssHeader.str().length();
  of << ssHeader.str();


  size_t n_cols = 0;
  size_t n_rows = 0;
  size_t n_entries = 0;
  std::vector<BUSData> v;
  v.reserve(N);
  uint64_t current_bc = 0xFFFFFFFFFFFFFFFFULL;
  //temporary data
  std::vector<int32_t> ecs;
  std::vector<int32_t> glist;
  ecs.reserve(100);
  std::vector<int32_t> u;
  u.reserve(100);
  std::vector<int32_t> column_v;
  std::vector<std::pair<int32_t, double>> column_vp;
  if (!opt.count_collapse) {
    column_v.reserve(N); 
  } else {
    column_vp.reserve(N);
    glist.reserve(100);
  }

  //Allocate histograms
  //Indexed as gene*histmax + histIndex
  size_t n_genes = genenames.size();
  const uint32_t histmax = 100;//set molecules with more than histmax copies to histmax 
  std::vector<double> histograms;
  std::vector<int> cellUMIs; //used for calculating CU per cell
  std::vector<int> cellCounts; //used for calculating CU per cell
  if (opt.count_gen_hist) {
	  
	  histograms = std::vector<double>(n_genes * histmax, 0);
	  cellUMIs.reserve(100000);
	  cellCounts.reserve(100000);
  }

  //set up random number generator for downsampling
  std::random_device						rand_dev;
  std::mt19937								generator(rand_dev());
  std::uniform_real_distribution<double>	distr(0.0, 1.0);


  //barcodes 
  std::vector<uint64_t> barcodes;
  int bad_count = 0;
  int compacted = 0;
  int rescued = 0;

  auto write_barcode_matrix = [&](const std::vector<BUSData> &v) {
    if(v.empty()) {
      return;
    }
    column_v.resize(0);
    n_rows+= 1;
    
    barcodes.push_back(v[0].barcode);
    double val = 0.0;
    size_t n = v.size();

    for (size_t i = 0; i < n; ) {
      size_t j = i+1;
      for (; j < n; j++) {
        if (v[i].UMI != v[j].UMI) {
          break;
        }
      }

      // v[i..j-1] share the same UMI
      ecs.resize(0);
      for (size_t k = i; k < j; k++) {
        ecs.push_back(v[k].ec);
      }
      
      if (opt.umi_gene_collapse) {
        intersect_genes_of_ecs(ecs,ec2genes, glist);
      }
      if (opt.umi_gene_collapse && glist.size() == 0) {
        // Gene-intersection zero, check for UMI collision
        std::vector<int32_t> ecs_within_molecule;
        while (ecs.size() > 0) {
          ecs_within_molecule.resize(0);
          for (size_t k = 0; k < ecs.size(); k++) {
            ecs_within_molecule.push_back(ecs[k]);
            intersect_genes_of_ecs(ecs_within_molecule, ec2genes, glist);
            if (glist.size() == 0) {
              ecs_within_molecule.pop_back();
            } else {
              ecs.erase(ecs.begin() + k);
              k--;
            }
          }
          int32_t ec = intersect_ecs(ecs_within_molecule, u, genemap, ecmap, ecmapinv, ec2genes);
          if (ec == -1) {
            ec = intersect_ecs_with_genes(ecs_within_molecule, genemap, ecmap, ecmapinv, ec2genes);
            if (ec == -1) {
              bad_count += ecs_within_molecule.size();
            } else {
              bool filter = false;
              if (!opt.count_gene_multimapping) {
                filter = (ec2genes[ec].size() != 1);
              }
              if (!filter) {
                rescued += ecs_within_molecule.size();
                column_v.push_back(ec);
              }
            }
          } else {
            bool filter = false;
            if (!opt.count_gene_multimapping) {
              filter = (ec2genes[ec].size() != 1);
            }
            if (!filter) {
              compacted += ecs_within_molecule.size()-1;
              column_v.push_back(ec);
            }
          }
        }
      } else { // Consider the records, which share the same UMI, as having come from the same molecule
        int32_t ec = intersect_ecs(ecs, u, genemap, ecmap, ecmapinv, ec2genes);
        if (ec == -1) {
          ec = intersect_ecs_with_genes(ecs, genemap, ecmap, ecmapinv, ec2genes);              
          if (ec == -1) {
            bad_count += j-i;
          } else {
            bool filter = false;
            if (!opt.count_gene_multimapping) {
              filter = (ec2genes[ec].size() != 1);
            }
            if (!filter) {
              rescued += j-i;
              column_v.push_back(ec);
            }
          }
          
        } else {
          bool filter = false;
          if (!opt.count_gene_multimapping) {
            filter = (ec2genes[ec].size() != 1);
          }
          if (!filter) {
            compacted += j-i-1;
            column_v.push_back(ec);
          }
        }
      }
      i = j; // increment
    }
    std::sort(column_v.begin(), column_v.end());
    size_t m = column_v.size();
    for (size_t i = 0; i < m; ) {
      size_t j = i+1;
      for (; j < m; j++) {
        if (column_v[i] != column_v[j]) {
          break;
        }
      }
      double val = j-i;
      of << n_rows << " " << (column_v[i]+1) << " " << val << "\n";
      n_entries++;
      
      i = j; // increment
    }
  };

  auto write_barcode_matrix_collapsed = [&](const std::vector<BUSData> &v) {
    if(v.empty()) {
      return;
    }
    column_vp.resize(0);
    n_rows+= 1;
	cellUMIs.push_back(0);
	cellCounts.push_back(0);

    barcodes.push_back(v[0].barcode);
    double val = 0.0;
    size_t n = v.size();

    std::vector<std::vector<int32_t>> ambiguous_genes;

    for (size_t i = 0; i < n; ) {
      size_t j = i+1;
      for (; j < n; j++) {
        if (v[i].UMI != v[j].UMI) {
          break;
        }
      }
	
      // v[i..j-1] share the same UMI
      ecs.resize(0);
	  uint32_t counts = 0;
      for (size_t k = i; k < j; k++) {
        ecs.push_back(v[k].ec);
		counts += v[k].count;
      }

      intersect_genes_of_ecs(ecs,ec2genes, glist);
      int gn = glist.size();
	  if (opt.count_downsampling_factor != 1.0) {
		  uint32_t newCounts = 0;
		  for (uint32_t c = 0; c < counts; ++c) {
			  if (distr(generator) <= opt.count_downsampling_factor) {
				  ++newCounts;
			  }
		  }
		  counts = newCounts;
		  if (newCounts == 0) {
			  gn = 0;//trick to skip quantification below
		  }

      }
      if (gn > 0) {
        if (opt.count_gene_multimapping) {
          for (auto x : glist) {
            column_vp.push_back({x, (opt.count_raw_counts ? counts : 1.0)/gn});
          }
		  
		  //Fill in histograms for prediction.
		  if (opt.count_gen_hist) {
			for (auto x : glist) {
				if (x < n_genes) { //crasches with an invalid gene file otherwise
					histograms[x * histmax + std::min(counts - 1, histmax - 1)] += 1.0 / gn; //histmax-1 since histograms[g][0] is the histogram value for 1 copy and so forth
				} else {
					std::cerr << "Mismatch between gene file and bus file, the bus file contains gene indices that is outside the gene range!\n";
				}
			}
			cellUMIs[barcodes.size()-1]++;
			cellCounts[barcodes.size()-1] += counts;
		  }
        } else {
          if (gn==1) {
            column_vp.push_back({glist[0],opt.count_raw_counts ? counts : 1.0});
			//Fill in histograms for prediction.
			if (opt.count_gen_hist) {
				if (glist[0] < n_genes) { //crasches with an invalid gene file otherwise
					histograms[glist[0] * histmax + std::min(counts - 1, histmax - 1)] += 1.0; //histmax-1 since histograms[g][0] is the histogram value for 1 copy and so forth
					cellUMIs[barcodes.size()-1]++;
					cellCounts[barcodes.size()-1] += counts;
				} else {
					std::cerr << "Mismatch between gene file and bus file, the bus file contains gene indices that is outside the gene range!\n";
				}
			}

          } else if (opt.count_em) {
            ambiguous_genes.push_back(std::move(glist));
			//Fill in histograms for prediction. This is a simplification. TODO: should be fixed for the em algorithm!
			if (opt.count_gen_hist) {
				for (auto x : glist) {
					if (x < n_genes) { //crasches with an invalid gene file otherwise
						histograms[x * histmax + std::min(counts - 1, histmax - 1)] += 1.0 / gn; //histmax-1 since histograms[g][0] is the histogram value for 1 copy and so forth
					} else {
						std::cerr << "Mismatch between gene file and bus file, the bus file contains gene indices that is outside the gene range!\n";
					}
				}
				cellUMIs[barcodes.size()-1]++;
				cellCounts[barcodes.size()-1] += counts;
			}
		  }
        }
      } else if (opt.umi_gene_collapse) { // Gene-intersection zero, check for UMI collision
        std::vector<int32_t> ecs_within_molecule;
        while (ecs.size() > 0) {
          ecs_within_molecule.resize(0);
          for (size_t k = 0; k < ecs.size(); k++) {
            ecs_within_molecule.push_back(ecs[k]);
            intersect_genes_of_ecs(ecs_within_molecule, ec2genes, glist);
            if (glist.size() == 0) {
              ecs_within_molecule.pop_back();
            } else {
              ecs.erase(ecs.begin() + k);
              k--;
            }
          }
          gn = glist.size();
          if (gn > 0) {
            if (opt.count_gene_multimapping) {
              for (auto x : glist) {
                column_vp.push_back({x, 1.0/gn});
              }
            } else {
              if (gn==1) {
                column_vp.push_back({glist[0],1.0});
              } else if (opt.count_em) {
                ambiguous_genes.push_back(std::move(glist));
              }
            }
          }
        }
      }
      i = j; // increment
    }
    std::sort(column_vp.begin(), column_vp.end());
    size_t m = column_vp.size();
    std::unordered_map<int32_t, double> col_map(m);
    std::vector<int32_t> cols;

    for (size_t i = 0; i < m; ) {
      size_t j = i+1;
      double val = column_vp[i].second;
      for (; j < m; j++) {
        if (column_vp[i].first != column_vp[j].first) {
          break;
        }
        val += column_vp[j].second;
      }
      col_map.insert({column_vp[i].first,val});
      cols.push_back(column_vp[i].first);

      n_entries++;
      
      i = j; // increment
    }

    if (opt.count_em) {
      //std::cerr << "Running EM algorithm" << std::endl;
      std::unordered_map<int32_t, double> c1,c2;
      // initialize with unique counts
      for (const auto &x : cols) {
        double val = 0;
        auto it = col_map.find(x);
        if (it != col_map.end()) {
          val = it->second;
        }
        c1.insert({x,val});
        c2.insert({x,0.0});
      }

      // c1 is from previous round, c2 is the new one
      int round = 0;
      double tol = 1e-3;
      double err = 1.0;
      if (!ambiguous_genes.empty()) {        
        while (round < 100 && err >= tol) {
          for (auto &x : col_map) {
            c2[x.first] = x.second;
          }

          for (const auto &glist : ambiguous_genes) {
            double s = 0.0;
            for (const auto &x : glist) {
              auto it = c1.find(x);
              if (it != c1.end()) {
                s += it->second;
              }
            }
            if (s > 0.0) {
              for (const auto &x : glist) {
                auto it = c1.find(x);
                if (it != c1.end()) {
                  c2[x] += it->second/s;
                }
              }
            }
          }

          err = 0.0;
          for (const auto &x : c1) {
            double a=x.second,b=c2[x.first];
            if (a>0 && b > 0) {
              err += std::abs((a-b)/(a+b));
            }
          }
          //std::cerr << "Round " << round << " error = " << err << std::endl;
          std::swap(c1,c2);
          round++;
        }
        std::swap(col_map,c1);
      }

    }



    for (const auto &x : cols) {
      double val = 0;
      auto it = col_map.find(x);
      if (it != col_map.end()) {
        val = it->second;
      }
      of << n_rows << " " << (x+1) << " " << val << "\n";
    }
  };

  for (const auto& infn : opt.files) { 
    std::streambuf *inbuf;
    std::ifstream inf;
    if (!opt.stream_in) {
      inf.open(infn.c_str(), std::ios::binary);
      inbuf = inf.rdbuf();
    } else {
      inbuf = std::cin.rdbuf();
    }
    std::istream in(inbuf); 

    parseHeader(in, h);
    bclen = h.bclen;
    
    int rc = 0;
    while (true) {
      in.read((char*)p, N*sizeof(BUSData));
      size_t rc = in.gcount() / sizeof(BUSData);
      nr += rc;
      if (rc == 0) {
        break;
      }

      
      for (size_t i = 0; i < rc; i++) {
        if (p[i].barcode != current_bc) {                 
          // output whatever is in v
          if (!v.empty()) {
            if (!opt.count_collapse) {
              write_barcode_matrix(v);
            } else {
				write_barcode_matrix_collapsed(v);
            }
          }
          v.clear();
          current_bc = p[i].barcode;
        }
        v.push_back(p[i]);

      }            
    }
    if (!v.empty()) {
      if (!opt.count_collapse) {
        write_barcode_matrix(v);
      } else {
        write_barcode_matrix_collapsed(v);
      }
    }

    if (!opt.stream_in) {
      inf.close();
    }
  }
  delete[] p; p = nullptr;

  if (!opt.count_collapse) {
    n_cols = ecmap.size();
  } else {
    n_cols = genenames.size();
  }

  of.close();
  
  //Rewrite header in a way that works for both Windows and Linux
  std::stringstream ss;
  ss << n_rows << " " << n_cols << " " << n_entries;
  std::string header = ss.str();
  int hlen = header.size();
  header = header + std::string(66 - hlen, ' ') + '\n';
  of.open(mtx_ofn, std::ios::in | std::ios::out);
  of << headerComments << header;
  of.close();

  // write updated ec file
  h.ecs = std::move(ecmap);
  if (!opt.count_collapse) {
    writeECs(ec_ofn, h);
  } else {
    writeGenes(gene_ofn, genenames);
  }
  // write barcode file
  std::ofstream bcof;
  bcof.open(barcodes_ofn);
  for (const auto &x : barcodes) {
    bcof << binaryToString(x, bclen) << "\n";
  }
  bcof.close();

  //write histogram file
  if (opt.count_gen_hist) {
	std::ofstream histof;
	histof.open(hist_ofn);

	for (size_t g = 0; g < genenames.size(); ++g) {
		//Indexed as gene*histmax + histIndex
		unsigned int offs = g * histmax;
		
		//first figure out the length of the histogram, don't write that to make the file smaller
		unsigned int histEnd = histmax - 1;
		for (; histEnd != 0; --histEnd) {
			if (histograms[offs + histEnd] != 0) {
				break;
			}
		}
		for (size_t c = 0; c <= histEnd; ++c) {
			if (c != 0) {
				histof << '\t';
			}
			histof << histograms[offs + c];
		}

		histof << "\n";
	}
	histof.close();
  }
  
  if (opt.count_gen_hist) {
	//write mean counts per UMI file (per gene)
	
	std::ofstream cuof;
	cuof.open(cu_ofn);
	//write header
	cuof << "gene\tCU\tUMIs\n"; 

	//prepare gene names for writing
	std::vector<std::string> names;
    names.resize(genenames.size());
	for (const auto &x : genenames) {
		if (x.second >= 0) {
			names[x.second] = x.first;
		}
	}


	for (size_t g = 0; g < genenames.size(); ++g) {
		//Indexed as gene*histmax + histIndex
		unsigned int offs = g * histmax;
		
		//calculate counts per UMI as the mean of the histogram
		double wsum = 0;
		double sum = 0;
		for (size_t c = 0; c < histmax; ++c) {
			wsum += double(c+1) * histograms[offs + c];
			sum += histograms[offs + c];
		}
		double cu = wsum/sum;
		if (sum == 0) {
			cuof << names[g] << '\t' << "NA" << '\t' << sum << '\n';
		} else {
			cuof << names[g] << '\t' << cu << '\t' << sum << '\n';
		}
	}
	cuof.close();
	
	//write cu per cell file
	
	std::ofstream cupcof;
	cupcof.open(cu_per_cell_ofn);
	//write header
	cupcof << "barcode\tCU\tUMIs\n"; 

	for (size_t bc = 0; bc < barcodes.size(); ++bc) {
		cupcof << binaryToString(barcodes[bc], bclen) << '\t' << double(cellCounts[bc]) / double(cellUMIs[bc]) << '\t' << cellUMIs[bc] << '\n';
	}
	cupcof.close();
  }
  

  //std::cerr << "bad counts = " << bad_count <<", rescued  =" << rescued << ", compacted = " << compacted << std::endl;

  //std::cerr << "Read in " << nr << " BUS records" << std::endl;
}

void bustools_count_mult(Bustools_opt &opt) {
  BUSHeader h;
  size_t nr = 0;
  size_t N = 100000;
  uint32_t bclen = 0;
  BUSData* p = new BUSData[N];

  // read and parse the equivalence class files

  std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> ecmapinv;
  std::vector<std::vector<int32_t>> ecmap;

  std::unordered_map<std::string, int32_t> txnames;
  parseTranscripts(opt.count_txp, txnames);
  std::vector<int32_t> genemap(txnames.size(), -1);
  std::unordered_map<std::string, int32_t> genenames;
  parseGenes(opt.count_genes, txnames, genemap, genenames);
  parseECs(opt.count_ecs, h);
  ecmap = std::move(h.ecs);
  ecmapinv.reserve(ecmap.size());
  for (int32_t ec = 0; ec < ecmap.size(); ec++) {
    ecmapinv.insert({ecmap[ec], ec});
  }
  std::vector<std::vector<int32_t>> ec2genes;        
  create_ec2genes(ecmap, genemap, ec2genes);


  std::ofstream of;
  std::string mtx_ofn = opt.output + ".mtx";
  std::string barcodes_ofn = opt.output + ".barcodes.txt";
  std::string ec_ofn = opt.output + ".ec.txt";
  std::string gene_ofn = opt.output + ".genes.txt";
  of.open(mtx_ofn); 

  // write out the initial header
  of << "%%MatrixMarket matrix coordinate real general\n%\n";
  // number of genes
  auto mat_header_pos = of.tellp();
  std::string dummy_header(66, '\n');
  for (int i = 0; i < 33; i++) {
    dummy_header[2*i] = '%';
  }
  of.write(dummy_header.c_str(), dummy_header.size());


  size_t n_cols = 0;
  size_t n_rows = 0;
  size_t n_entries = 0;
  std::vector<BUSData> v;
  v.reserve(N);
  uint64_t current_bc = 0xFFFFFFFFFFFFFFFFULL;
  //temporary data
  std::vector<int32_t> ecs;
  std::vector<int32_t> glist;
  ecs.reserve(100);
  std::vector<int32_t> u;
  u.reserve(100);
  std::vector<int32_t> column_v;
  std::vector<std::pair<int32_t, double>> column_vp;
  if (!opt.count_collapse) {
    column_vp.reserve(N); 
  } else {
    column_vp.reserve(N);
    glist.reserve(100);
  }
  //barcodes 
  std::vector<uint64_t> barcodes;
  int bad_count = 0;
  int compacted = 0;
  int rescued = 0;


  auto write_barcode_matrix = [&](const std::vector<BUSData> &v) {
    if(v.empty()) {
      return;
    }
    column_vp.resize(0);
    n_rows+= 1;
    
    barcodes.push_back(v[0].barcode);
    double val = 0.0;
    size_t n = v.size();
    
    for (size_t i = 0; i < n; i++) {
      int32_t ec = v[i].ec;
      if (!opt.count_gene_multimapping) {
        ecs.resize(0);
        ecs.push_back(ec);
        intersect_genes_of_ecs(ecs, ec2genes, glist);
        int gn = glist.size();
        if (gn != 1) {
          continue;
        }
      }
      column_vp.push_back({ec,v[i].count});
    }
    std::sort(column_vp.begin(), column_vp.end());
    size_t m = column_vp.size();
    for (size_t i = 0; i < m; ) {
      size_t j = i+1;
      double val = column_vp[i].second;
      for (; j < m; j++) {
        if (column_vp[i].first != column_vp[j].first) {
          break;
        }
        val += column_vp[j].second;
      }
      n_entries++;
      of << n_rows << " " << (column_vp[i].first+1) << " " << val << "\n";
      i = j; // increment
    }
  };
  
  auto write_barcode_matrix_collapsed = [&](const std::vector<BUSData> &v) {
    if(v.empty()) {
      return;
    }
    column_vp.resize(0);
    n_rows+= 1;
    
    barcodes.push_back(v[0].barcode);
    double val = 0.0;
    size_t n = v.size();

    for (size_t i = 0; i < n; i++) {
      ecs.resize(0);
      ecs.push_back(v[i].ec);
      
      intersect_genes_of_ecs(ecs, ec2genes, glist);
      int gn = glist.size();
      if (gn > 0) {
        if (opt.count_gene_multimapping) {
          for (auto x : glist) {
            column_vp.push_back({x, v[i].count/gn});
          }
        } else {
          if (gn==1) {
            column_vp.push_back({glist[0],v[i].count});
          } 
        }
      }
    }

    std::sort(column_vp.begin(), column_vp.end());
    size_t m = column_vp.size();
    std::unordered_map<int32_t, double> col_map(m);
    std::vector<int32_t> cols;

    for (size_t i = 0; i < m; ) {
      size_t j = i+1;
      double val = column_vp[i].second;
      for (; j < m; j++) {
        if (column_vp[i].first != column_vp[j].first) {
          break;
        }
        val += column_vp[j].second;
      }
      col_map.insert({column_vp[i].first,val});
      cols.push_back(column_vp[i].first);

      n_entries++;
      
      i = j; // increment
    }



    for (const auto &x : cols) {
      double val = 0;
      auto it = col_map.find(x);
      if (it != col_map.end()) {
        val = it->second;
      }
      of << n_rows << " " << (x+1) << " " << val << "\n";
    }

  };

  for (const auto& infn : opt.files) { 
    std::streambuf *inbuf;
    std::ifstream inf;
    if (!opt.stream_in) {
      inf.open(infn.c_str(), std::ios::binary);
      inbuf = inf.rdbuf();
    } else {
      inbuf = std::cin.rdbuf();
    }
    std::istream in(inbuf); 

    parseHeader(in, h);
    bclen = h.bclen;
    
    int rc = 0;
    while (true) {
      in.read((char*)p, N*sizeof(BUSData));
      size_t rc = in.gcount() / sizeof(BUSData);
      nr += rc;
      if (rc == 0) {
        break;
      }

      
      for (size_t i = 0; i < rc; i++) {
        if (p[i].barcode != current_bc) {                 
          // output whatever is in v
          if (!v.empty()) {
            if (!opt.count_collapse) {
              write_barcode_matrix(v);
              } else {
                write_barcode_matrix_collapsed(v);
              }
            }
          v.clear();
          current_bc = p[i].barcode;
        }
        v.push_back(p[i]);

      }            
    }
    if (!v.empty()) {
      if (!opt.count_collapse) {
        write_barcode_matrix(v);
      } else {
        write_barcode_matrix_collapsed(v);
      }
    }

    if (!opt.stream_in) {
      inf.close();
    }
  }
  delete[] p; p = nullptr;

  if (!opt.count_collapse) {
    n_cols = ecmap.size();
  } else {
    n_cols = genenames.size();
  }

  of.close();
  
  std::stringstream ss;
  ss << n_rows << " " << n_cols << " " << n_entries << "\n";
  std::string header = ss.str();
  int hlen = header.size();
  assert(hlen < 66);
  of.open(mtx_ofn, std::ios::binary | std::ios::in | std::ios::out);
  of.seekp(mat_header_pos);
  of.write("%",1);
  of.write(std::string(66-hlen-2,' ').c_str(),66-hlen-2);
  of.write("\n",1);
  of.write(header.c_str(), hlen);
  of.close();

  // write updated ec file
  h.ecs = std::move(ecmap);
  if (!opt.count_collapse) {
    writeECs(ec_ofn, h);
  } else {
    writeGenes(gene_ofn, genenames);
  }
  // write barcode file
  std::ofstream bcof;
  bcof.open(barcodes_ofn);
  for (const auto &x : barcodes) {
    bcof << binaryToString(x, bclen) << "\n";
  }
  bcof.close();
  //std::cerr << "bad counts = " << bad_count <<", rescued  =" << rescued << ", compacted = " << compacted << std::endl;

  //std::cerr << "Read in " << nr << " BUS records" << std::endl;
}