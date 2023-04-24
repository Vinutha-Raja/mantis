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

#include "sparsepp/spp.h"
#include "tsl/sparse_map.h"

#include "MantisFS.h"
#include "ProgOpts.h"
#include "coloreddbg.h"
#include "squeakrconfig.h"
#include "json.hpp"
#include "mantis_utils.hpp"
#include "mantisconfig.hpp"

/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:  
 * ============================================================================
 */
	int
build_main ( BuildOpts& opt )
{
	spdlog::logger* console = opt.console.get();
    std::cout << opt.inlist;
	std::ifstream infile(opt.inlist);
  uint64_t num_samples{0};
  if (infile.is_open()) {
    std::string line;
    while (std::getline(infile, line)) { ++num_samples; }
    infile.clear();
    infile.seekg(0, std::ios::beg);
    console->info("Will build mantis index over {} input experiments.", num_samples);
  } else {
    console->error("Input file {} does not exist or could not be opened.", opt.inlist);
    std::exit(1);
  }

  /** try and create the output directory
   *  and write a file to it.  Complain to the user
   *  and exit if we cannot.
   **/
	std::string prefix(opt.out);
	if (prefix.back() != '/') {
		prefix += '/';
	}
	// make the output directory if it doesn't exist
	if (!mantis::fs::DirExists(prefix.c_str())) {
		mantis::fs::MakeDir(prefix.c_str());
	}
	// check to see if the output dir exists now
	if (!mantis::fs::DirExists(prefix.c_str())) {
		console->error("Output dir {} could not be successfully created.", prefix);
		exit(1);
	}

  std::vector<SampleObject<CQF<KeyObject>*>> inobjects;
  std::vector<CQF<KeyObject>> cqfs;

	// reserve QF structs for input CQFs
  inobjects.reserve(num_samples);
  cqfs.reserve(num_samples);

    // lmer vs list of file names containing the lmer
    std::unordered_map<std::string, std::list<int>> lmerMap;

	// mmap all the input cqfs
	std::string squeakr_file;
    std::unordered_map<std::string, int> fileIDMapping;
    int fileID = 1;
    while (infile >> squeakr_file) {
        std::string line;
        fileIDMapping[squeakr_file] = fileID;
        std::ifstream sqfile(squeakr_file);
        while (std::getline(sqfile, line)) {
            lmerMap[line].push_back(fileIDMapping[squeakr_file]);
        }
        fileID++;
    }

	uint32_t nqf = 0;
	uint32_t kmer_size{0};
	console->info("Reading input Squeakr files.");

  {
    std::ofstream jfile(prefix + "/" + mantis::LMER_INDEX_FILE);
    std::ofstream sfile(prefix + "/" + mantis::SAMPLEID_FILE);
    if (jfile.is_open()) {
        for (const auto& entry : lmerMap) {
            jfile << entry.first << ":";
            for (const auto& value : entry.second) {
                jfile << " " << value;
            }
            jfile << "\n";
        }

    } else {
      console->error("Could not write to output directory {}", prefix);
    }
    jfile.close();

    if (sfile.is_open()) {
          for (const auto& entry : fileIDMapping) {
              sfile << entry.first << ":" << entry.second<<std::endl;
          }

      } else {
          console->error("Could not write to sample ID file {}", prefix);
      }
      sfile.close();

  }

  return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
