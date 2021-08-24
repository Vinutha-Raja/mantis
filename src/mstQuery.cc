//
// Created by Fatemeh Almodaresi on 2018-10-15.
//
#include <fstream>
#include <vector>
#include <CLI/Timer.hpp>
#include <canonicalKmer.h>

#include "ProgOpts.h"
#include "kmer.h"
#include "mstQuery.h"


void MSTQuery::loadMST(std::string indexDir) {
    sdsl::load_from_file(parentbv, indexDir + mantis::PARENTBV_FILE);
    sdsl::load_from_file(deltabv, indexDir + mantis::DELTABV_FILE);
    sdsl::load_from_file(bbv, indexDir + mantis::BOUNDARYBV_FILE);
    sbbv = sdsl::bit_vector::select_1_type(&bbv);
    zero = parentbv.size() - 1; // maximum color id which
    logger->info("Loaded the new color class index");
    logger->info("\t--> parent size: {}", parentbv.size());
    logger->info("\t--> delta size: {}", deltabv.size());
    logger->info("\t--> boundary size: {}", bbv.size());
}

std::vector<uint64_t> MSTQuery::buildColor(uint64_t eqid, QueryStats &queryStats,
                                           LRUCacheMap *lru_cache,
                                           RankScores *rs,
                                           std::unordered_map<uint64_t, std::vector<uint64_t>> *fixed_cache,
                                           nonstd::optional<uint64_t> &toDecode) {

    (void) rs;
    std::vector<uint32_t> flips(numSamples);
    std::vector<uint32_t> xorflips(numSamples, 0);
    uint64_t i{eqid}, from{0}, to{0};
    int64_t height{0}, weight{0};
    auto froms = std::vector<uint64_t>();
    queryStats.totEqcls++;
    bool foundCache = false;
    uint32_t iparent = parentbv[i];
    while (iparent != i) {
        if (fixed_cache and fixed_cache->find(i) != fixed_cache->end()) {
            const auto &vs = (*fixed_cache)[i];
            for (auto v : vs) {
                xorflips[v] = 1;
            }
            queryStats.cacheCntr++;
            foundCache = true;
            break;
        } else if (lru_cache and lru_cache->contains(i)) {
            const auto &vs = (*lru_cache)[i];
            for (auto v : vs) {
                xorflips[v] = 1;
            }
            queryStats.cacheCntr++;
            foundCache = true;
            break;
        }
        from = (i > 0) ? (sbbv(i) + 1) : 0;
//        std::cerr << " -> " << i;
        froms.push_back(from);

        if (queryStats.trySample) {
            if ((!toDecode) and
//                (occ > 10) and
                (height > 10)
                and
                (lru_cache and
                 !lru_cache->contains(iparent))) {
                        toDecode = iparent;
            }
        }
        i = iparent;
        iparent = parentbv[i];
        ++queryStats.totSel;
        ++height;
    }
    if (!foundCache and i != zero) {
//        std::cerr << " -> " << i;
        from = (i > 0) ? (sbbv(i) + 1) : 0;
        froms.push_back(from);
        ++queryStats.totSel;
        queryStats.rootedNonZero++;
        ++height;
    }
//    std::cerr << "\n";
    for (auto f : froms) {
        queryStats.noCacheCntr++;
        bool found = false;
        uint64_t wrd{0};
        auto start = f;
//        std::cerr << " -> " << start << " : ";
        do {
            uint64_t wrdLen = std::min(static_cast<uint64_t >(64), bbv.size()-start);
            wrd = bbv.get_int(start, wrdLen);
            for (uint64_t j = 0; j < wrdLen; j++) {
//                std::cerr << deltabv[start+j] << " ";
                flips[deltabv[start + j]] ^= 0x01;
                weight++;
                if ((wrd >> j) & 0x01) {
                    found = true;
                    break;
                }
            }
            start += wrdLen;
        } while (!found);
//        std::cerr << "\n";
    }

    /*if (foundCache) {
        queryStats.heightDist.push_back(height);
        queryStats.weightDist.push_back(weight);
    } else {
        queryStats.noCache_heightDist.push_back(height);
        queryStats.noCache_weightDist.push_back(weight);
    }*/
    std::vector<uint64_t> eq;
    eq.reserve(numWrds);
    for (i = 0; i < numSamples; i++) {
        if (flips[i] ^ xorflips[i]) {
            eq.push_back(i);
        }
    }
    return eq;
}

void MSTQuery::findSamples(ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> &cdbg,
                           LRUCacheMap &lru_cache,
                           RankScores *rs,
                           QueryStats &queryStats,
                           uint64_t numQueries) {
    uint64_t ksize{cdbg.get_current_cqf()->keybits()}, numBlocks{cdbg.get_numBlocks()};
    std::vector<std::unordered_set<mantis::KmerHash>> blockKmers(numBlocks);
    std::vector<std::ofstream> blockKmerFiles;
    std::vector<uint64_t> blockKmerCnt(numBlocks, 0);
    for (uint64_t i = 0; i < numBlocks; i++) {
        blockKmerFiles.emplace_back(outfile+".tmp" + std::to_string(i), std::ios::out | std::ios::binary);
    }
    // split kmers based on minimizers into blocks
    logger->info("# of distinct kmers to query: {}", kmer2cidMap.size());
    for (auto kv : kmer2cidMap) {
        auto minimizers = cdbg.findMinimizer(kv.first, ksize); //assuming not hashed
        blockKmers[cdbg.minimizerBlock[minimizers.first]].insert(kv.first);
        blockKmerFiles[cdbg.minimizerBlock[minimizers.first]].write(reinterpret_cast<const char*>(&kv.first), sizeof(kv.first));
        blockKmerCnt[cdbg.minimizerBlock[minimizers.first]]++;
    }

    for (auto &f : blockKmerFiles) f.close();

    // go block by block and query kmers
    std::vector<uint64_t> query_colors;
    query_colors.reserve(10000000);
    for (uint64_t blockId = 0; blockId < numBlocks; blockId++) {
        if (blockKmerCnt[blockId] == 0)
            continue;
        logger->info("{}: kmer count: {}", blockId, blockKmers[blockId].size());
        cdbg.replaceCQFInMemory(blockId);
        logger->info("\t-->cqf loaded");
        std::vector<uint64_t> kmers(blockKmerCnt[blockId]);
        std::ifstream inputf(outfile+".tmp" + std::to_string(blockId), std::ios::in | std::ios::binary);
        inputf.read(reinterpret_cast<char *>(kmers.data()), sizeof(uint64_t)*blockKmerCnt[blockId]);
        inputf.close();
        std::string cmd = "rm " + outfile+".tmp" + std::to_string(blockId);
        system(cmd.c_str());
        for (auto key : kmers) {
            KeyObject keyObj(key, 0, 0);
            uint64_t eqclass = cdbg.query_kmerInCurDbg(keyObj, 0);
            if (eqclass) {
                kmer2cidMap[key].first = eqclass - 1;
                query_colors.push_back(eqclass - 1);
            } else {
                kmer2cidMap[key].first = invalid;
            }
        }
        std::sort(query_colors.begin(), query_colors.end());
        // tbb::parallel_sort(query_colors.begin(), query_colors.end());
        query_colors.erase(std::unique(query_colors.begin(), query_colors.end()), query_colors.end());
    }


    cdbg.replaceCQFInMemory(invalid);
    auto colorIt = boomphf::range(query_colors.begin(), query_colors.end());
    colorMph = boomphf::mphf<uint64_t, boomphf::SingleHashFunctor<uint64_t>>(query_colors.size(), colorIt, threadCount, 2.5);

    std::vector<uint64_t> eqQueriesStartIdx(query_colors.size()+1,0), eqQueriesEndIdx(query_colors.size()+1,0);
    std::vector<uint32_t> eqQueries;
    for (auto &kv : kmer2cidMap) {
        if (kv.second.first != invalid) {
            eqQueriesEndIdx[colorMph.lookup(kv.second.first)] += kv.second.second.size();
        }
    }
    for (uint64_t i = 1; i < eqQueriesEndIdx.size(); i++) {
        eqQueriesEndIdx[i] += eqQueriesEndIdx[i-1];
        eqQueriesStartIdx[i] = eqQueriesEndIdx[i-1];
    }
    eqQueries.resize(eqQueriesEndIdx.back());
    for (auto &kv : kmer2cidMap) {
        if (kv.second.first != invalid) {
            auto idx = colorMph.lookup(kv.second.first);
            for (auto &q : kv.second.second) {
                eqQueries[eqQueriesStartIdx[idx]] = q;
                eqQueriesStartIdx[idx]++;
            }
        }
    }
    eqQueriesStartIdx[0] = 0;
    for (uint64_t i = 1; i < eqQueriesEndIdx.size(); i++) {
        if (eqQueriesStartIdx[i] != eqQueriesEndIdx[i]) {
            std::cerr << "Error in calculating the start and end for queries containing a color\n";
            std::cerr << eqQueriesStartIdx[i] << "!=" << eqQueriesEndIdx[i] << "\n";
            std::exit(3);
        }
        eqQueriesStartIdx[i] = eqQueriesEndIdx[i-1];
    }
    // we can now clear kmer2cidMap

    nonstd::optional<uint64_t> toDecode{nonstd::nullopt};
    nonstd::optional<uint64_t> dummy{nonstd::nullopt};
    logger->info("Fetching colors now ..");
    if (!MSTIsLoaded) {
        loadMST(prefix);
    }
    allQueries.resize(numQueries);
    for (auto &r : allQueries) {
        r.resize(numSamples+1, 0); // last index keeps count of distinct kmers
    }
    for (auto &it : query_colors) {
        uint64_t eqclass_id = it;

        std::vector<uint64_t> setbits;
        if (lru_cache.contains(eqclass_id)) {
            setbits = lru_cache[eqclass_id];//.get(eqclass_id);
//            setbits = (*lru_cache[eqclass_id]);//.get(eqclass_id);
            queryStats.cacheCntr++;
//            queryStats.heightDist.push_back(0);
//            queryStats.weightDist.push_back(0);
        } else {
            queryStats.noCacheCntr++;
            toDecode.reset();
            dummy.reset();
            queryStats.trySample = (queryStats.noCacheCntr % 10 == 0);
            setbits = buildColor(eqclass_id, queryStats, &lru_cache, rs, nullptr, toDecode);
            lru_cache.emplace(eqclass_id, setbits);
            if ((queryStats.trySample) and toDecode) {
                auto s = buildColor(*toDecode, queryStats, nullptr, nullptr, nullptr, dummy);
                lru_cache.emplace(*toDecode, s);
            }
        }

        // update query sample counts
        auto eqMappedIndx = colorMph.lookup(eqclass_id);
        for (auto idx = eqQueriesStartIdx[eqMappedIndx]; idx < eqQueriesEndIdx[eqMappedIndx]; idx++) {
            auto q = eqQueries[idx];
            allQueries[q][numSamples]++;
            for (auto &c : setbits) {
                allQueries[q][c]++;
            }
        }
    }
}

uint64_t MSTQuery::parseBulkKmers(const std::string &file, u_int64_t kmer_size) {
    uint64_t numOfQueries{0};
    std::string read;
    std::ifstream ipfile(file);
    while (ipfile >> read) {
        parseKmers(numOfQueries, read, indexK);
        numOfQueries++;
    }
    ipfile.close();
    for (auto &kv : kmer2cidMap) {
        auto &v = kv.second.second;
        std::sort(v.begin(), v.end());
        // tbb::parallel_sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    return numOfQueries;
}

void MSTQuery::parseKmers(uint32_t readId, std::string read, uint64_t kmer_size) {
    //CLI::AutoTimer timer{"First round going over the file ", CLI::Timer::Big};
//    std::cerr << "\r" << readId << " " << read.length() << "   ";
    bool done = false;
    std::vector<uint32_t> emptyVec;
    while (!done and read.length() >= kmer_size) {
        uint64_t first = 0;
        uint64_t first_rev = 0;
        uint64_t item = 0;
        bool allNuclValid = true;
        for (uint32_t i = 0; i < kmer_size; i++) { //First kmer
            uint8_t curr = Kmer::map_base(read[i]);
            if (curr > DNA_MAP::G) { // 'N' is encountered
                if (i + 1 < read.length())
                    read = read.substr(i + 1);
                else
                    done = true; // you've reached the end of a sequence (last nucl. is N)
                allNuclValid = false;
                break;
            }
            first = first | curr;
            first = first << 2;
        }
        if (allNuclValid) {
            first = first >> 2;
            first_rev = static_cast<uint64_t >(Kmer::reverse_complement(first, kmer_size));

            //cout << "kmer: "; cout << int_to_str(first);
            //cout << " reverse-comp: "; cout << int_to_str(first_rev) << endl;

            if (Kmer::compare_kmers(first, first_rev))
                item = first;
            else
                item = first_rev;
            if (kmer2cidMap.find(item) == kmer2cidMap.end()) {
                kmer2cidMap[item] = std::make_pair(std::numeric_limits<uint64_t>::max(), emptyVec);
            }
            kmer2cidMap[item].second.push_back(readId);

            uint64_t next = (first << 2) & BITMASK(2 * kmer_size);
            uint64_t next_rev = first_rev >> 2;
            uint64_t i = 0;
            for (i = kmer_size; i < read.length(); i++) { //next kmers
                //cout << "K: " << read.substr(i-K+1,K) << endl;
                uint8_t curr = Kmer::map_base(read[i]);
                if (curr > DNA_MAP::G) { // 'N' is encountered
                    if (i + 1 < read.length())
                        read = read.substr(i + 1);
                    else
                        done = true;
                    break;
                }
                next |= curr;
                auto tmp = static_cast<uint64_t>(Kmer::reverse_complement_base(curr));
                tmp <<= (kmer_size * 2 - 2);
                next_rev = next_rev | tmp;
                if (Kmer::compare_kmers(next, next_rev))
                    item = next;
                else
                    item = next_rev;
                if (kmer2cidMap.find(item) == kmer2cidMap.end()) {
                    kmer2cidMap[item] = std::make_pair(std::numeric_limits<uint64_t>::max(), emptyVec);
                }
                kmer2cidMap[item].second.push_back(readId);
                next = (next << 2) & BITMASK(2 * kmer_size);
                next_rev = next_rev >> 2;
            }
            if (i == read.length()) done = true;
        }
    }
}

void MSTQuery::reset() {
    kmer2cidMap.clear();
}

void MSTQuery::clear() {
    reset();
    parentbv.resize(0);
    deltabv.resize(0);
    bbv.resize(0);
}

void output_results(MSTQuery &mstQuery,
                    std::ofstream &opfile,
                    const std::vector<std::string> &sampleNames,
                    QueryStats &queryStats,
                    uint64_t numQueries) {
    //CLI::AutoTimer timer{"Second round going over the file + query time ", CLI::Timer::Big};
    mantis::QueryResults &result = mstQuery.allQueries;//mstQuery.getResultList(numQueries);
    for (auto &q : result) {
        // last element in the result for each query contains # of distinct kmers
        opfile << "seq" << queryStats.cnt++ << '\t' << q[q.size()-1] << '\n';
        for (uint64_t i = 0; i < q.size()-1; i++) {
            if (q[i] > 0) {
                opfile << sampleNames[i] << '\t' << q[i] << '\n';
            }
        }
    }
}

void output_results_json(MSTQuery &mstQuery,
                         std::ofstream &opfile,
                         std::vector<std::string> &sampleNames,
                         QueryStats &queryStats,
                         uint64_t nquery,
                         bool isFirstQuery) {
    //CLI::AutoTimer timer{"Query time ", CLI::Timer::Big};
//    mantis::QueryResults result = mstQuery.getResultList(nquery);
    mantis::QueryResults &result = mstQuery.allQueries;//mstQuery.getResultList(numQueries);
    for (auto &q : result) {
        if (!isFirstQuery) {
            opfile << ",";
            isFirstQuery = false;
        }
        opfile << "{ \"qnum\": " << queryStats.cnt++ << ",  \"num_kmers\": "
               << q[q.size() - 1] << ", \"res\": {\n";
        bool isFirst = true;
        for (uint64_t i = 0; i < q.size()-1; i++) {
            if (q[i] > 0) {
                if (!isFirst) {
                    opfile << ",\n";
                }
                opfile << " \"" << sampleNames[i] << "\": " << q[i];
                isFirst = false;
            }
        }
        opfile << "}}\n";
    }
}

std::vector<std::string> loadSampleFile(const std::string &sampleFileAddr) {
    std::vector<std::string> sampleNames;
    std::ifstream sampleFile(sampleFileAddr);
    uint64_t id;
    std::string name;
    while (sampleFile >> id >> name) {
        if (id >= sampleNames.size()) {
            sampleNames.resize(id + 1);
        }
        sampleNames[id] = name;
    }
    return sampleNames;
}

/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:
 * ============================================================================
 */
int mst_query_main(QueryOpts &opt) {
    uint32_t queryK = opt.k;
    QueryStats queryStats;

    spdlog::logger *logger = opt.console.get();
    std::string sample_file(opt.prefix + mantis::SAMPLEID_FILE);

    std::vector<std::string> sampleNames = loadSampleFile(sample_file);
    queryStats.numSamples = sampleNames.size();
    logger->info("Number of experiments: {}", queryStats.numSamples);

    logger->info("Loading the first cdbg and the first CQF...");
    ColoredDbg<SampleObject<CQF<KeyObject> *>, KeyObject> cdbg(opt.prefix, MANTIS_DBG_IN_MEMORY);
    cdbg.set_console(logger);
    auto curDbg = cdbg.get_current_cqf();
    if (not curDbg) {
        logger->error("No dbg has been loaded into memory.");
        std::exit(3);
    }
    uint64_t indexK = curDbg->keybits() / 2;
    if (queryK == 0) queryK = indexK;
    logger->info("Done loading cqf. k is {}", indexK);

    logger->info("Loading color classes in the MST form...");
    MSTQuery mstQuery(opt.prefix, opt.output, indexK, queryK, queryStats.numSamples, logger, !opt.process_in_bulk);
    logger->info("Done Loading color classes. Total # of color classes is {}",
                 mstQuery.parentbv.size() - 1);

    logger->info("Querying colored dbg.");
    std::ofstream opfile(opt.output);
    LRUCacheMap cache_lru(100000);
    RankScores rs(1);
    uint64_t numOfQueries{0};
    CLI::AutoTimer timer{"query time ", CLI::Timer::Big};
    if (opt.process_in_bulk) {
        numOfQueries = mstQuery.parseBulkKmers(opt.query_file, indexK);
        logger->info("# of observed sequences: {}", numOfQueries);
        mstQuery.findSamples(cdbg, cache_lru, &rs, queryStats, numOfQueries);
        if (opt.use_json) {
            opfile << "[\n";
            output_results_json(mstQuery, opfile, sampleNames, queryStats, numOfQueries, true);
            opfile << "]\n";
        } else {
            output_results(mstQuery, opfile, sampleNames, queryStats, numOfQueries);
        }
    } else {
        std::ifstream ipfile(opt.query_file);
        std::string read;
        if (opt.use_json) {
            bool isFirstQuery = true;
            opfile << "[\n";
            while (ipfile >> read) {
                mstQuery.reset();
                mstQuery.parseKmers(numOfQueries, read, indexK);
                mstQuery.findSamples(cdbg, cache_lru, &rs, queryStats, 1);
                output_results_json(mstQuery, opfile, sampleNames, queryStats, 1, isFirstQuery);
                isFirstQuery = false;
                numOfQueries++;
            }
            opfile << "]\n";
        } else {
            while (ipfile >> read) {
                mstQuery.reset();
                mstQuery.parseKmers(numOfQueries, read, indexK);
                mstQuery.findSamples(cdbg, cache_lru, &rs, queryStats, 1);
                output_results(mstQuery, opfile, sampleNames, queryStats, 1);
                numOfQueries++;
            }
        }
        ipfile.close();
    }
    opfile.close();
    logger->info("Writing done.");

    logger->info("cache was used {} times and not used {} times",
                 queryStats.cacheCntr, queryStats.noCacheCntr);
    logger->info("total selects = {}, time per select = {}",
                 queryStats.totSel, queryStats.selectTime.count() / queryStats.totSel);
    logger->info("total # of queries = {}, total # of queries rooted at a non-zero node = {}",
                 queryStats.totEqcls, queryStats.rootedNonZero);
    return EXIT_SUCCESS;
}