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

	std::string dbg_file(prefix + mantis::CQF_FILE);
	std::string sample_file(prefix + mantis::SAMPLEID_FILE);
	std::vector<std::string> eqclass_files = mantis::fs::GetFilesExt(prefix.c_str(),
                                                                   mantis::EQCLASS_FILE);

	ColoredDbg<SampleObject<CQF<KeyObject>*>, KeyObject> cdbg(dbg_file,
																														eqclass_files,
																														sample_file,
																														MANTIS_DBG_IN_MEMORY);
	uint64_t kmer_size = cdbg.get_cqf()->keybits() / 2;
  console->info(" 1 Read colored dbg with {} k-mers and {} color classes",
                cdbg.get_cqf()->dist_elts(), cdbg.get_num_bitvectors());

	//cdbg.get_cqf()->dump_metadata(); 
	//CQF<KeyObject> cqf(query_file, false);
	//CQF<KeyObject>::Iterator it = cqf.begin(1);
	//mantis::QuerySet input_kmers;
	//do {
		//KeyObject k = *it;
		//input_kmers.insert(k.key);
		//++it;
	//} while (!it.done());

	//mantis::QuerySets multi_kmers;
	//multi_kmers.push_back(input_kmers);

	console->info("Reading query kmers from disk.");
	uint32_t seed = 2038074743;
	uint64_t total_kmers = 0;
    std::unordered_map<mantis::KmerHash, uint64_t> uniqueKmers;
	mantis::QuerySets multi_kmers = Kmer::parse_kmers(query_file.c_str(),
																										kmer_size,
																										total_kmers,
																										opt.process_in_bulk,
																										uniqueKmers);
	console->info("Total k-mers to query: {}", total_kmers);
    console->info("Number of uniqueKmers to query: {}", uniqueKmers.size());
    for (auto const &pair: uniqueKmers) {
        console->info("Kmer: {}, Count: {}", pair.first, pair.second);
    }
    console->info("Number of multi_kmers to query: {}", multi_kmers.size());

    //Convert Kmers to Lmers
    uint64_t lmer_size = 15;
    MinimizerScanner scanner(kmer_size, lmer_size, 0, true, 0);
    mantis::QuerySets multi_lmers;
    std::unordered_map<mantis::KmerHash, uint64_t> uniqueLmers;
    for (auto const &pair: uniqueKmers) {
        string seq = Kmer::int_to_str(pair.first, kmer_size);
        scanner.LoadSequence(seq);
        uint64_t *mmp = scanner.NextMinimizer();
        uniqueLmers[*mmp] = pair.second;
//        console->info("kmer: {}, count: {}", pair.first, pair.second);
    }

    for(auto element : multi_kmers){
        std::unordered_set<uint64_t> lmerSet;
        for (auto it = element.begin(); it != element.end(); ++it) {
            string seq = Kmer::int_to_str(*it, kmer_size);
            cout << "kmer: "; cout << Kmer::int_to_str(*it, kmer_size);
            scanner.LoadSequence(seq);
            uint64_t *mmp = scanner.NextMinimizer();
            cout<< "lmer: "<<*mmp<<" lmer str:"<<Kmer::int_to_str(*mmp, lmer_size)<<std::endl;
            lmerSet.insert(*mmp);
        }
        multi_lmers.push_back(lmerSet);
    }

	std::ofstream opfile(output_file);
	console->info("Querying the colored dbg.");

  if (use_json) {
      // Passing lmers instead of kmers
    output_results_json(multi_lmers, cdbg, opfile, opt.process_in_bulk, uniqueLmers);
  } else {
      // Passing lmers instead of kmers
    output_results(multi_lmers, cdbg, opfile, opt.process_in_bulk, uniqueLmers);
  }
	//std::cout << "Writing samples and abundances out." << std::endl;
	opfile.close();
	console->info("Writing done.");

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
