/*
 * ============================================================================
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * ============================================================================
 */


#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <unordered_set>
#include <set>
#include <bitset>
#include <cassert>
#include <fstream>
#include <list>

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <openssl/rand.h>

#include "MantisFS.h"
#include "ProgOpts.h"
#include "spdlog/spdlog.h"
#include "kmer.h"
#include "coloreddbg.h"
#include "common_types.h"
#include "CLI/CLI.hpp"
#include "CLI/Timer.hpp"
#include "mantisconfig.hpp"
#include "minimizer.h"

void output_results(mantis::QuerySets& multi_kmers,
										ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject>&
										cdbg, std::ofstream& opfile, bool is_bulk,
                    std::unordered_map<mantis::KmerHash, uint64_t> &uniqueKmers) {
  mantis::QueryResults qres;
	uint32_t cnt= 0;
  {
    CLI::AutoTimer timer{"Query time ", CLI::Timer::Big};
    if (is_bulk) {
        std::unordered_map<uint64_t, std::vector<uint64_t>> result = cdbg.find_samples(uniqueKmers);
        for (auto& kmers : multi_kmers) {
            opfile <<  cnt++ << '\t' << kmers.size() << '\n';
            std::vector<uint64_t> kmerCnt(cdbg.get_num_samples());
            for (auto &k : kmers) {
                if (uniqueKmers[k]) {
                    for (auto experimentId : result[k]) {
                        kmerCnt[experimentId]++;
                    }
                }
            }
            for (auto i=0; i<kmerCnt.size(); ++i) {
                if (kmerCnt[i] > 0)
                    opfile << cdbg.get_sample(i) << '\t' << kmerCnt[i] << '\n';
            }
            //++qctr;
        }
    } else {
        for (auto &kmers : multi_kmers) {
            //std::sort(kmers.begin(), kmers.end());
            opfile << cnt++ << '\t' << kmers.size() << '\n';
            mantis::QueryResult result = cdbg.find_samples(kmers);
            for (auto i = 0; i < result.size(); ++i) {
                if (result[i] > 0)
                    opfile << cdbg.get_sample(i) << '\t' << result[i] << '\n';
            }
            //++qctr;
        }
    }
  }
}

void output_lmer_results(mantis::QuerySets& multi_kmers, std::string samplefile, std::string dbgfile, std::ofstream& opfile, bool is_bulk,
                    std::unordered_map<mantis::KmerHash, uint64_t> &uniqueKmers) {
    mantis::QueryResults qres;
    std::ifstream infile(samplefile);
    std::ifstream indexfile(dbgfile);

    // Reading File ID mapping function
    std::unordered_map<int, std::string> fileIDMapping;
    std::string line;
    while (std::getline(infile, line)) {
        size_t pos = line.find(':');
        std::string value = line.substr(0, pos);
        int key = std::stoi(line.substr(pos + 1));
        fileIDMapping[key] = value;
    }

    // Reading the mantis index file function
    std::unordered_map<uint64_t , std::list<int>> lmerIndex;

    while (std::getline(indexfile, line)) {
        std::stringstream ss(line);
        std::string key;
        std::getline(ss, key, ':');
        uint64_t lmer = Kmer::str_to_int(key);
//        cout<< "map"<< key << " " << lmer <<std::endl;

        int id;
        while (ss >> id) {
            lmerIndex[lmer].push_back(id);
        }
    }

    uint32_t cnt= 0;
    {
        CLI::AutoTimer timer{"Query time ", CLI::Timer::Big};
        if (is_bulk) {

            for (auto& kmers : multi_kmers) {
                opfile <<  cnt++ << '\t' << kmers.size() << '\n';
                std::unordered_map<int, int> kmerCnt;
                for (auto &k : kmers) {
                        for (auto experimentId : lmerIndex[k]) {
                            kmerCnt[experimentId]++;
                        }
                }
                for (auto i=0; i<kmerCnt.size(); ++i) {
                    if (kmerCnt[i] > 0)
                        opfile << fileIDMapping[i] << '\t' << kmerCnt[i] << '\n';
                }
            }
        } else {
            for (auto &kmers : multi_kmers) {
                opfile << cnt++ << '\t' << kmers.size() << '\n';
                std::unordered_map<int, int> kmerCnt;
                for (auto &k : kmers) {
//                    cout << k << " " << Kmer::int_to_str(k, 7) << endl;
                    for (auto experimentId : lmerIndex[k]) {
                        kmerCnt[experimentId]++;
                    }
                }
                for (auto i=0; i<kmerCnt.size(); ++i) {
                    if (kmerCnt[i] > 0)
                        opfile << fileIDMapping[i] << '\t' << kmerCnt[i] << '\n';
                }
            }
        }
    }
}

void output_results_json(mantis::QuerySets& multi_kmers,
												 ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject>&
												 cdbg, std::ofstream& opfile, bool is_bulk,
                         std::unordered_map<mantis::KmerHash, uint64_t> &uniqueKmers) {
  mantis::QueryResults qres;
	uint32_t cnt= 0;
  {
    CLI::AutoTimer timer{"Query time ", CLI::Timer::Big};

    opfile << "[\n";
      size_t qctr{0};
      size_t nquery{multi_kmers.size()};

      if (is_bulk) {
          std::unordered_map<uint64_t, std::vector<uint64_t>> result = cdbg.find_samples(uniqueKmers);
          for (auto& kmers : multi_kmers) {
              opfile << "{ \"qnum\": " << cnt++ << ",  \"num_kmers\": " << kmers.size() << ", \"res\": {\n";
              std::vector<uint64_t> kmerCnt(cdbg.get_num_samples());
              for (auto &k : kmers) {
                  if (uniqueKmers[k]) {
                      for (auto experimentId : result[k]) {
                          kmerCnt[experimentId]++;
                      }
                  }
              }
              for (auto i=0; i<kmerCnt.size(); ++i) {
                  if (kmerCnt[i] > 0) {
                          opfile << " \"" << cdbg.get_sample(i) << "\": " << kmerCnt[i];
                      if (i != kmerCnt.size()) {
                          opfile << ",\n";
                      }
                  }
              }
              opfile << "}}";
              if (qctr < nquery - 1) { opfile << ","; }
              opfile << "\n";
              ++qctr;
          }
      }
      else {
          for (auto &kmers : multi_kmers) {
              //std::sort(kmers.begin(), kmers.end());
              opfile << "{ \"qnum\": " << cnt++ << ",  \"num_kmers\": " << kmers.size() << ", \"res\": {\n";
              mantis::QueryResult result = cdbg.find_samples(kmers);
              uint64_t sampleCntr = 0;
              for (auto it = result.begin(); it != result.end(); ++it) {
                  if (*it > 0)
                      opfile << " \"" << cdbg.get_sample(sampleCntr) << "\": " << *it;
                  if (std::next(it) != result.end()) {
                      opfile << ",\n";
                  }
                  sampleCntr++;
              }
              opfile << "}}";
              if (qctr < nquery - 1) { opfile << ","; }
              opfile << "\n";
              ++qctr;
          }
      }
    opfile << "]\n";
  }
}


/* 
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:  
 * ============================================================================
 */
int query_main (QueryOpts& opt)
{
  //CLI::App app("Mantis query");

  std::string prefix = opt.prefix;
  uint64_t lmer_size = opt.l;
  uint64_t kmer_size = opt.k;
  std::string query_file = opt.query_file;
  std::string output_file = opt.output;//{"samples.output"};
  bool use_json = opt.use_json;

  // Make sure the prefix is a full folder
  if (prefix.back() != '/') {
    prefix.push_back('/');
  }
  // make the output directory if it doesn't exist
  if (!mantis::fs::DirExists(prefix.c_str())) {
    mantis::fs::MakeDir(prefix.c_str());
  }

  spdlog::logger* console = opt.console.get();
	console->info("Reading colored dbg from disk.");

	std::string dbg_file(prefix + mantis::LMER_INDEX_FILE);
	std::string sample_file(prefix + mantis::SAMPLEID_FILE);

	console->info("Reading query kmers from disk.");
	uint32_t seed = 2038074743;
	uint64_t total_kmers = 0;
    std::unordered_map<mantis::KmerHash, uint64_t> uniqueKmers;

    //Convert Kmers to Lmers
//    MinimizerScanner scanner(kmer_size, lmer_size, 0, true, 0);
	mantis::QuerySets multi_kmers = Kmer::parse_kmers(query_file.c_str(),
																										kmer_size,
																										total_kmers,
																										opt.process_in_bulk,
																										uniqueKmers, lmer_size);
	console->info("Total k-mers to query: {}", total_kmers);
    
	std::ofstream opfile(output_file);
	console->info("Querying the colored dbg.");

  if (use_json) {
      // Passing lmers instead of kmers
//    output_results_json(multi_kmers, cdbg, opfile, opt.process_in_bulk, uniqueKmers);
  } else {
      // Passing lmers instead of kmers
      output_lmer_results(multi_kmers, sample_file, dbg_file,  opfile, opt.process_in_bulk, uniqueKmers);
  }
	//std::cout << "Writing samples and abundances out." << std::endl;
	opfile.close();
	console->info("Writing done.");

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
