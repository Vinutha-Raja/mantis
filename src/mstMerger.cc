#include <memory>

#include <utility>

//
// Created by Fatemeh Almodaresi.
//
#include <string>
#include <sstream>
#include <cstdio>
#include <stdio.h>
#include <unordered_map>
#include <iterator>
#include <algorithm>
#include <zconf.h>
#include <cqfMerger.h>

#include "MantisFS.h"
#include "mstMerger.h"
#include "ProgOpts.h"
#include "mstPlanner.h"


#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) - 1ULL)

static constexpr uint64_t MAX_ALLOWED_TMP_EDGES_IN_FILE{16000000};//{130000000}; // Up to 2G
static constexpr uint64_t SHIFTBITS{32};
static constexpr uint64_t HIGHBIT_MASK{(1ULL << 32) - 1ULL};
static constexpr uint64_t LOWBIT_MASK{std::numeric_limits<uint64_t>::max() - HIGHBIT_MASK};
static constexpr uint64_t CONSTNUMBlocks{1 << 16};

MSTMerger::MSTMerger(std::string prefixIn, spdlog::logger *loggerIn, uint32_t numThreads,
                     std::string prefixIn1, std::string prefixIn2) :
                prefix(std::move(prefixIn)),
                nThreads(numThreads), numBlocks(CONSTNUMBlocks) {
    logger = loggerIn;//.get();
    prefixes[0] = prefixIn1;
    prefixes[1] = prefixIn2;
    // Make sure the prefix is a full folder
    if (prefix.back() != '/') {
        prefix.push_back('/');
    }
    if (prefixes[0].back() != '/') {
        prefixes[0].push_back('/');
    }
    if (!mantis::fs::DirExists(prefixes[0].c_str())) {
        logger->error("Index parent directory for first mst, {}, does not exist", prefixes[0]);
        std::exit(1);
    }
    if (prefixes[1].back() != '/') {
        prefixes[1].push_back('/');
    }
    if (!mantis::fs::DirExists(prefixes[1].c_str())) {
        logger->error("Index parent directory for second mst, {}, does not exist", prefixes[1]);
        std::exit(1);
    }

    for (auto &q:queryStats)
        q.resize(nThreads);

    for (auto &lru: lru_cache) {
        for (uint64_t t = 0; t < nThreads; t++) {
            lru.emplace_back(1000);
        }
    }

    std::string sample_file = prefixes[0] + mantis::SAMPLEID_FILE;
    std::ifstream sampleid(sample_file);
    std::string tmp;
    while (sampleid >> tmp >> tmp) {
        numSamples++;
    }
    sampleid.close();
    toBeMergedNumOfSamples[0] = numSamples; // this line is important!

    sample_file = prefixes[1] + mantis::SAMPLEID_FILE;//(prefix.c_str() , mantis::SAMPLEID_FILE);
    sampleid.open(sample_file);
    while (sampleid >> tmp >> tmp) {
        numSamples++;
    }
    toBeMergedNumOfSamples[1] = numSamples - toBeMergedNumOfSamples[0];
    sampleid.close();
    maxWeightInFile = std::min(static_cast<uint32_t >(maxWeightInFile), numSamples);

    uint64_t newColorIdCnt, cIdx, c1, c2;
    std::ifstream cp(prefix + "newID2oldIDs");
    // colorPairs colors for mantis 1 and mantis 2 are 0-based --> color ID starts from 0
    cp.read(reinterpret_cast<char *>(&newColorIdCnt), sizeof(newColorIdCnt));
    cp.read(reinterpret_cast<char *>(&ccCnt[0]), sizeof(ccCnt[0]));
    cp.read(reinterpret_cast<char *>(&ccCnt[1]), sizeof(ccCnt[1]));
    auto c1len{static_cast<uint64_t >(ceil(log2(ccCnt[0])))}, c2len{static_cast<uint64_t >(ceil(log2(ccCnt[1])))};
    logger->info("# of color classes based on count of colorPairs: {}", newColorIdCnt);
    colorPairs[0] = sdsl::int_vector<>(newColorIdCnt, 0, c1len);
    colorPairs[1] = sdsl::int_vector<>(newColorIdCnt, 0, c2len);
    uint64_t maxIndex = 0;
    for (auto i = 0; i < newColorIdCnt; i++) {
        cp.read(reinterpret_cast<char *>(&cIdx), sizeof(cIdx));
        cp.read(reinterpret_cast<char *>(&c1), sizeof(c1));
        cp.read(reinterpret_cast<char *>(&c2), sizeof(c2));
        colorPairs[0][cIdx] = c1;
        colorPairs[1][cIdx] = c2;
        maxIndex = maxIndex>=cIdx?maxIndex:cIdx;
    }
    cp.close();
    num_colorClasses = colorPairs[0].size() + 1;
    zero = num_colorClasses-1;

    logger->info("# of experiments: {}", numSamples);
}

/**
 * Builds an MST consists of 3 main steps:
 * 1. construct the color graph for all the colorIds derived from dbg
 *      This phase just requires loading the CQF
 * 2. calculate the weights of edges in the color graph
 *      This phase requires at most two buffers of color classes
 * 3. find MST of the weighted color graph
 */
void MSTMerger::mergeMSTs() {
    auto t_start = time(nullptr);
    buildEdgeSets();
    calculateMSTBasedWeights();
    encodeColorClassUsingMST();

//    std::string cmd = "rm " + prefix + "newID2oldIDs";
//    system(cmd.c_str());
    auto t_end = time(nullptr);
    logger->info("MST merge completed in {} s.", t_end - t_start);
}

/**
 * iterates over all elements of CQF,
 * find all the existing neighbors, and build a color graph based on that
 * @return pair of maximum observed colorId and total observed edges (including duplicates)
 */
bool MSTMerger::buildEdgeSets() {

    std::vector<uint64_t> cnts(nThreads, 0);
    uint64_t maxId{0}, numOfKmers{0};
    std::vector<std::string> cqfBlocks = mantis::fs::GetFilesExt(prefix.c_str(), mantis::CQF_FILE);
    std::sort(cqfBlocks.begin(), cqfBlocks.end());
    std::vector<std::pair<uint64_t, uint64_t>> tmpEdges;
    tmpEdges.reserve(MAX_ALLOWED_TMP_EDGES_IN_FILE);

    uint64_t sampledEdges = 10000;
    spp::sparse_hash_map<std::pair<uint64_t , uint64_t >, uint64_t, Custom_Pair_Hasher> popularEdges;
    /*std::string cqf_file(cqfBlocks[0]);
    CQF<KeyObject> cqf(cqf_file, CQF_FREAD);
    k = cqf.keybits() / 2;
    auto it = cqf.begin();
    while (!it.done() and popularEdges.size() < sampledEdges) {
        KeyObject keyObject = *it;
        std::vector<std::pair<uint64_t, uint64_t>> edgeList;
        spp::sparse_hash_map<std::pair<uint64_t , uint64_t >, uint64_t, Custom_Pair_Hasher> tmpPopEdg;
        findNeighborEdges(cqf, keyObject, edgeList, tmpPopEdg);
        for (auto &v : edgeList) {
            if (popularEdges.find(v) == popularEdges.end()) {
                popularEdges[v] = 0;
            }
            popularEdges[v]++;
        }
        ++it;
    }
    std::cerr << "Before doubling popular edges\n";
    usleep(10000000);

    logger->info("PopularEdges");
    for (auto &kv : popularEdges) {
        tmpEdges.push_back(kv.first);
    }
    std::cerr << "After doubling popular edges\n";
    usleep(10000000);*/
    for (uint64_t c = 0; c < cqfBlocks.size(); c++) {
        logger->info("Reading colored dbg from disk...");
//        std::cerr << cqfBlocks[c] << "\n";
        std::cerr << "\n\ncqf" << c << "\n";
        std::string cqf_file(cqfBlocks[c]);
        CQF<KeyObject> cqf(cqf_file, CQF_FREAD);
        k = cqf.keybits() / 2;
        cqf.dump_metadata();
        std::cerr << "\n\n";
        logger->info("Done loading cdbg. k is {}, numThreads is {}", k, nThreads);
        logger->info("Iterating over cqf & building edgeSet ...");
        // build color class edges in a multi-threaded manner
        std::vector<std::thread> threads;
        for (uint32_t i = 0; i < nThreads; ++i) {
            threads.emplace_back(std::thread(&MSTMerger::buildPairedColorIdEdgesInParallel, this, i,
                                             std::ref(cqf), std::ref(tmpEdges), std::ref(curFileIdx),
                                             std::ref(cnts[i]),
                                             std::ref(maxId), std::ref(numOfKmers),
                                             std::ref(popularEdges)));
        }
        for (auto &t : threads) { t.join(); }
    }
    writeMutex.lock();
    if (not tmpEdges.empty()) {
        std::cerr << "tmpEdges.size(): " << tmpEdges.size() << " sizeof(element): " << sizeof(std::remove_reference<decltype(tmpEdges)>::type::value_type) <<  "\n";
        edgePairSortUniq(tmpEdges);
        std::string filename(prefix+"tmp"+std::to_string(curFileIdx));
        std::ofstream tmpfile(filename, std::ios::out | std::ios::binary);
        uint64_t maxIdd = 0;
        std::cerr << "writing to file tmp" << curFileIdx << "\n";
        std::cerr << "Total elements: " << tmpEdges.size() << " size:"
                  << sizeof(std::remove_reference<decltype(tmpEdges)>::type::value_type)*tmpEdges.size()
                  << " with maxcolorID: " << maxIdd << "\n";
        uint64_t vecSize = tmpEdges.size();
        totalWrittenEdges = vecSize;
        tmpfile.write(reinterpret_cast<const char *>(&vecSize), sizeof(vecSize));
        tmpfile.write(reinterpret_cast<const char *>(tmpEdges.data()),
                      sizeof(std::remove_reference<decltype(tmpEdges)>::type::value_type)*tmpEdges.size());
        tmpfile.close();
        curFileIdx++;
    }
    writeMutex.unlock();
    logger->info("Done writing the last tmp file");
    for (uint32_t i = 0; i < nThreads; ++i) {
        totalWrittenEdges += cnts[i];
    }
    num_colorClasses = maxId + 1;
    zero = static_cast<colorIdType>(num_colorClasses);
    num_colorClasses++;// last one is for zero as the dummy color class with ID equal to actual num of color classes
    logger->info("Total number of kmers and color classes (including dummy) observed: {}, {}.\nTotal number of edges (including duplicates): {}", numOfKmers, num_colorClasses, totalWrittenEdges);
    return true;
}

void MSTMerger::buildPairedColorIdEdgesInParallel(uint32_t threadId,
                                            CQF<KeyObject> &cqf,
                                            std::vector<std::pair<uint64_t, uint64_t>> &tmpEdges,
                                            uint64_t &curFileIdx,
                                            std::uint64_t& cnt,
                                            uint64_t &maxId, uint64_t &numOfKmers,
                                                  spp::sparse_hash_map<std::pair<uint64_t , uint64_t >, uint64_t, Custom_Pair_Hasher > &popularEdges) {
    //std::cout << "THREAD ..... " << threadId << " " << cqf.range() << "\n";
    uint64_t kmerCntr{0}, localMaxId{0};
    __uint128_t startPoint = threadId * (cqf.range() / (__uint128_t) nThreads);
    __uint128_t endPoint =
            threadId + 1 == nThreads ? cqf.range() + 1 : (threadId + 1) * (cqf.range() / (__uint128_t) nThreads);
    /*std::cerr << threadId << ": s" << (uint64_t) (startPoint/(__uint128_t)0xFFFFFFFFFFFFFFFF) << " "
              << "sr" << (uint64_t) (startPoint%(__uint128_t)0xFFFFFFFFFFFFFFFF) << " "
            << "e" << (uint64_t) (endPoint/(__uint128_t)0xFFFFFFFFFFFFFFFF) << " "
            << "er" << (uint64_t) (endPoint%(__uint128_t)0xFFFFFFFFFFFFFFFF) << "\n";*/
    auto tmpEdgeListSize = MAX_ALLOWED_TMP_EDGES_IN_FILE / nThreads;
//    auto tmpEdgeListSize = 10000;
    std::vector<std::pair<uint64_t , uint64_t >> edgeList;
    edgeList.reserve(tmpEdgeListSize);
    auto appendStore = [&]() {
        edgePairSortUniq(edgeList);
        colorMutex.lock();
        tmpEdges.insert(tmpEdges.end(), edgeList.begin(), edgeList.end());
        edgeList.clear();
        if (tmpEdges.size() >= MAX_ALLOWED_TMP_EDGES_IN_FILE) {
            std::vector<std::pair<uint64_t , uint64_t >> toWrite = std::move(tmpEdges);
            tmpEdges.clear();
            tmpEdges.shrink_to_fit();
            writeMutex.lock();
            colorMutex.unlock();
            edgePairSortUniq(toWrite);
            std::string filename(prefix + "tmp" + std::to_string(curFileIdx));
            std::ofstream tmpfile(filename, std::ios::out | std::ios::binary);
            uint64_t toWriteSize = toWrite.size();
            cnt += toWriteSize;
            tmpfile.write(reinterpret_cast<const char *>(&toWriteSize), sizeof(toWriteSize));
            tmpfile.write(reinterpret_cast<const char *>(toWrite.data()),
                          sizeof(std::remove_reference<decltype(toWrite)>::type::value_type)*toWrite.size());
            tmpfile.close();
            curFileIdx++;
            writeMutex.unlock();
        } else {
            colorMutex.unlock();
        }
    };
    auto it = cqf.setIteratorLimits(startPoint, endPoint);
    while (!it.reachedHashLimit()) {
        KeyObject keyObject = *it;
        uint64_t curEqId = keyObject.count - 1;
        //nodes[curEqId] = 1; // set the seen color class id bit
        localMaxId = curEqId > localMaxId ? curEqId : localMaxId;
        // Add an edge between the color class and each of its neighbors' colors in dbg
        findNeighborEdges(cqf, keyObject, edgeList, popularEdges);
        if (edgeList.size() >= tmpEdgeListSize) {
            appendStore();
        }
        ++it;
        kmerCntr++;
        /*if (kmerCntr % 10000000 == 0) {
            std::cerr << "\rthread " << threadId << ": Observed " << (numOfKmers + kmerCntr) / 1000000 << "M kmers and " << cnt << " edges";
        }*/
    }
    appendStore();
    colorMutex.lock();
    maxId = localMaxId > maxId ? localMaxId : maxId;
    numOfKmers += kmerCntr;
    std::cerr << "\r";
    logger->info("Thread {}: Observed {} kmers and {} edges", threadId, numOfKmers, cnt);
    colorMutex.unlock();
}

/**
 * loads the color class table in parts
 * calculate the hamming distance between the color bitvectors fetched from color class table
 * for each pair of color IDs
 * having w buckets where w is the maximum possible weight (number of experiments)
 * put the pair in its corresponding bucket based on the hamming distance value (weight)
 * @return true if successful
 */
bool MSTMerger::calculateMSTBasedWeights() {

    uint64_t mstZero[2];
    mstZero[0] = ccCnt[0] - 1;
    mstZero[1] = ccCnt[1] - 1;
    auto c1len{static_cast<uint64_t >(ceil(log2(ccCnt[0])))}, c2len{static_cast<uint64_t >(ceil(log2(ccCnt[1])))};

    std::cerr << "Running the MST planner\n";
//    usleep(10000000);
    mst[0] = std::make_unique<MSTQuery>(prefixes[0], k, k, toBeMergedNumOfSamples[0], logger);
    mst[1] = std::make_unique<MSTQuery>(prefixes[1], k, k, toBeMergedNumOfSamples[1], logger);


    auto mstPlanner = std::make_unique<MSTPlanner>(prefix,
                          mst[0].get(), mst[1].get(),
                          colorPairs[0], colorPairs[1],
                          curFileIdx,
                          MAX_ALLOWED_TMP_EDGES_IN_FILE/(2.0*curFileIdx),
                          nThreads,
                          logger,
                          true);

    for (auto mstIdx = 0; mstIdx < 2; mstIdx++) {
        mstPlanner->plan(mstIdx, mst[mstIdx].get(), fixed_cache[mstIdx]);
    }
    mstPlanner.reset(nullptr);


    logger->info("Pop counting for the two Manti color class vectors of size {}, {} respectively", ccCnt[0], ccCnt[1]);
    sdsl::int_vector<> ccSetBitCnts[2];
    for (auto mstIdx = 0; mstIdx < 2; ++mstIdx) {
        ccSetBitCnts[mstIdx] = sdsl::int_vector<>(ccCnt[mstIdx], 0, ceil(log2(toBeMergedNumOfSamples[mstIdx])));
        for (auto partition = 0; partition < 2; ++partition) {
            std::vector<std::thread> threads;
            for (auto t = 0; t < nThreads; ++t) {
                threads.emplace_back([&, mstIdx, t, partition]() {
                    uint64_t s = ((2*t+partition)/(double)(2*nThreads))*ccSetBitCnts[mstIdx].size();
                    uint64_t e = ((2*t+partition+1)/(double)(2*nThreads))*ccSetBitCnts[mstIdx].size();
                     for (auto i = s; i < e; i++) {
                         std::vector<uint64_t> eq;
                         buildMSTBasedColor(i, mst[mstIdx].get(), lru_cache[mstIdx][t], eq, queryStats[mstIdx][t],
                                            fixed_cache[mstIdx]);
                         ccSetBitCnts[mstIdx][i] = eq.size();
                     }
                 });
            }
            for (auto &t : threads) {
                t.join();
            }
            threads.clear();
           /* std::cerr << "mstIdx: " << mstIdx <<
                      " cacheCntr: " << queryStats[mstIdx][0].cacheCntr <<
                      " noCacheCntr: " << queryStats[mstIdx][0].noCacheCntr << "\n";*/
        }
    }

    logger->info("Opening {} output weight bucket files", maxWeightInFile);
    std::vector<uint64_t> bufferEndIndices;
    uint64_t size = splitHarmonically(MAX_ALLOWED_TMP_EDGES_IN_FILE/2, maxWeightInFile, bufferEndIndices);
    std::vector<std::pair<uint64_t , uint64_t >> output_buffer(size);
    std::vector<uint64_t> bufferCurrIndices(bufferEndIndices.size(), 0);
    std::vector<uint64_t > bufferCnt(bufferEndIndices.size(), 0);
    for (auto i = 1; i < bufferCurrIndices.size(); i++) {
        bufferCurrIndices[i] = bufferEndIndices[i-1];
    }
    std::vector<std::ofstream> outWeightFile;
    outWeightFile.reserve(maxWeightInFile);
    for (auto tmpOutIdx = 1; tmpOutIdx <= maxWeightInFile; tmpOutIdx++) {
        std::string filename = prefix+"w"+std::to_string(tmpOutIdx);
        outWeightFile.emplace_back(filename, std::ios::out | std::ios::binary);
        uint64_t cnt = 0;
        outWeightFile.back().write(reinterpret_cast<char*>(&cnt), sizeof(cnt));
    }

    //    std::vector<std::tuple<uint64_t , uint64_t , uint32_t >> inMemEdgeWeight;
    //    inMemEdgeWeight.reserve(1000000);
    // a lambda function to store the edges in the corresponding weight bucket file
    // managing its own output buffer
    auto store_edge_weight = [&](auto pair, auto w) {
        if (w <= maxWeightInFile) {
            w--;
            if (bufferCurrIndices[w] == bufferEndIndices[w]) {
                auto startIdx = w == 0?0:bufferEndIndices[w-1];
                outWeightFile[w].write(reinterpret_cast<char *>(output_buffer.data()+startIdx),
                                       sizeof(std::remove_reference<decltype(output_buffer)>::type::value_type) * (bufferEndIndices[w]-startIdx));
                bufferCnt[w] += (bufferEndIndices[w]-startIdx);
                bufferCurrIndices[w] = startIdx;
            }
            output_buffer[bufferCurrIndices[w]] = pair;
            bufferCurrIndices[w]++;
        }
        /*else {
            inMemEdgeWeight.emplace_back(pair.first, pair.second, w);
        }*/
    };

    logger->info("Loading {} temp edge files.", curFileIdx);
    uint64_t tmpedgecnt=0;
    std::vector<TmpFileIterator> tis;
    Minheap_edge minheap;
    for (uint64_t i = 0; i < curFileIdx; i++) {
        tis.emplace_back(prefix+"tmp"+std::to_string(i), MAX_ALLOWED_TMP_EDGES_IN_FILE/(2.0*curFileIdx));
    }
    for (uint64_t i = 0; i < curFileIdx; i++) {
        if (tis[i].end()) continue;
        minheap.push(&tis[i]);
    }
    uint64_t maxAllowedCnt = nThreads * 1000000;
    std::vector<std::pair<colorIdType , colorIdType>> uniqueEdges[2];
    std::vector<std::pair<colorIdType , colorIdType >> outputMstEdges;
    outputMstEdges.reserve(maxAllowedCnt);
    uniqueEdges[0].reserve(maxAllowedCnt);
    uniqueEdges[1].reserve(maxAllowedCnt);

    logger->info("Calculating weights of non-dummy edges... ");
//    usleep(10000000);
    auto val = minheap.top()->get_val();
    while (!minheap.empty()) {
        do {
            auto cur = minheap.top();
            val = cur->get_val();
            if (cur->next()) {
                minheap.replace_top(cur);
            }
            else
                minheap.pop();
            tmpedgecnt++;
        } while (!minheap.empty() && val == minheap.top()->get_val());
        outputMstEdges.push_back(val);
        // calculate weight
        for (uint64_t mstIdx = 0; mstIdx < 2; mstIdx++) {
            auto edge = std::make_pair(colorPairs[mstIdx][val.first], colorPairs[mstIdx][val.second]);
            if (edge.first != edge.second)
                uniqueEdges[mstIdx].push_back(edge);
        }
        if (outputMstEdges.size() > maxAllowedCnt or minheap.empty()) {
            std::unique_ptr<boomphf::mphf<std::pair<colorIdType , colorIdType>, Custom_Pair_Hasher>> mph[2];
            std::vector<uint32_t> weights[2];
            // call parallel edge calculator for both sides
            for (uint64_t mstIdx = 0; mstIdx < 2; mstIdx++) {
                edgePairSortUniq(uniqueEdges[mstIdx]);
                auto er = boomphf::range(uniqueEdges[mstIdx].begin(), uniqueEdges[mstIdx].end());
                mph[mstIdx] = std::make_unique<boomphf::mphf<std::pair<colorIdType , colorIdType>,
                        Custom_Pair_Hasher>>(uniqueEdges[mstIdx].size(), er, nThreads, 2.5, false);
               /* logger->info("Total memory consumed by all the MPH tables = {} MB.",
                              (mph[mstIdx]->totalBitSize() / 8) / (1024 * 1024));*/
                weights[mstIdx].resize(uniqueEdges[mstIdx].size(), 0);
                std::vector<std::thread> threads;
                for (auto threadId = 0; threadId < nThreads; ++threadId) {
                    threads.emplace_back([&, mstIdx, threadId]() {
                        uint64_t startIdx = (threadId/(double)nThreads)*uniqueEdges[mstIdx].size(),
                        endIdx = ((threadId+1)/(double)nThreads)*uniqueEdges[mstIdx].size();
                        for (auto idx = startIdx; idx < endIdx; ++idx) {
                            auto edge = uniqueEdges[mstIdx][idx];
                            if (edge.first == mstZero[mstIdx]) {
                                weights[mstIdx][mph[mstIdx]->lookup(edge)] = ccSetBitCnts[mstIdx][edge.second];
                            } else if (edge.second == mstZero[mstIdx]) {
                                weights[mstIdx][mph[mstIdx]->lookup(edge)] = ccSetBitCnts[mstIdx][edge.first];
                            } else {
                                weights[mstIdx][mph[mstIdx]->lookup(edge)] = mstBasedHammingDist(edge.first,
                                                                                                 edge.second,
                                                                                                 mst[mstIdx].get(),
                                                                                                 lru_cache[mstIdx][threadId],
                                                                                                 queryStats[mstIdx][threadId],
                                                                                                 fixed_cache[mstIdx]);
                            }
                        }
                    });
                }
                for (auto &t : threads) {
                    t.join();
                }
            }
            // go over the edges and store them based on weight:
            for (auto &e : outputMstEdges) {
                uint64_t ws[2] = {0,0};
                for (auto mstIdx = 0; mstIdx < 2; mstIdx++) {
                    auto edge = std::make_pair(colorPairs[mstIdx][e.first], colorPairs[mstIdx][e.second]);
                    if (edge.first != edge.second) {
                        ws[mstIdx] = weights[mstIdx][mph[mstIdx]->lookup(edge)];
                    }
                }
                auto w = ws[0]+ws[1];
                if (w == 0) {
                    logger->error("Weight cannot be zero for edge {}:({},{}) -> {}:({},{})",
                            e.first, colorPairs[0][e.first], colorPairs[1][e.first],
                            e.second, colorPairs[0][e.second], colorPairs[1][e.second]);
                    std::exit(3);
                }
                store_edge_weight(e, w);
            }
            logger->info("all: {}, uniq1: {}, uniq2:{}", outputMstEdges.size(), uniqueEdges[0].size(), uniqueEdges[1].size());
            outputMstEdges.clear();
            uniqueEdges[0].clear();
            uniqueEdges[1].clear();
        }
    }
    outputMstEdges.shrink_to_fit();
    uniqueEdges[0].shrink_to_fit();
    uniqueEdges[1].shrink_to_fit();

    // calculate and store weights of the edges connected to dummy node
    logger->info("Calculating weights of dummy edges... ");
    ccBitsBucketCnt.resize(numSamples, 0);
    if (maxWeightInFile < ccBitsBucketCnt.size()) {
        for (auto i = 0; i < colorPairs[0].size(); i++) {
            auto cc1 = colorPairs[0][i];
            auto cc2 = colorPairs[1][i];
            if (cc1 > mstZero[0] or cc2 > mstZero[1]) {
                logger->error("Should not happen. Either of the colorIDs does not exist. cc1={}, max1={}. cc2={}, max2={}",
                              cc1, mstZero[0], cc2, mstZero[1]);
            }
            auto w = ccSetBitCnts[0][cc1] + ccSetBitCnts[1][cc2];
            if (w == 0) {
                logger->error("The weight of an edge is zero! edge:{},{}->dummy", cc1,cc2);
                std::exit(3);
            }
            ccBitsBucketCnt[w-1]++;
        }
        ccBitsBucketCnt[maxWeightInFile-1] = 0;
        for (uint64_t i = maxWeightInFile; i < ccBitsBucketCnt.size(); i++) {
            ccBitsBucketCnt[i] += ccBitsBucketCnt[i - 1];
//            std::cerr << i << ":" << ccBitsBucketCnt[i] << "\n";
        }
    }
//    usleep(10000000);
    ccBits = sdsl::int_vector<>(ccBitsBucketCnt.back(), 0, ceil(log2(zero)));
    std::cerr << "After initializing the ccBits vector.\n";
//    usleep(10000000);
    for (auto i = 0; i < colorPairs[0].size(); i++) {
        auto cc1 = colorPairs[0][i];
        auto cc2 = colorPairs[1][i];
        auto w = ccSetBitCnts[0][cc1] + ccSetBitCnts[1][cc2];
        if (w == 0) {
            logger->error("The weight of an edge is zero! edge:{}->dummy or <{},{}>->dummy", i, cc1,cc2);
            std::exit(3);
        }
        if (w <= maxWeightInFile) {
            auto pair = std::make_pair(i, zero);
            store_edge_weight(pair, w);
        } else {
            auto &idx = ccBitsBucketCnt[w-1];
            idx--;
            ccBits[idx] = i;
        }
    }
    logger->info("Done calculating weight of all the edges!");

    // store remaining buffer to file for each weight
    uint64_t totalUniqEdgeCnt = 0;
    for (auto w = 0; w < outWeightFile.size(); w++) {
        auto startIdx = w == 0?0:bufferEndIndices[w-1];
        if (bufferCurrIndices[w] != startIdx) {
            outWeightFile[w].write(reinterpret_cast<char *>(output_buffer.data()+startIdx),
                                   sizeof(std::remove_reference<decltype(output_buffer)>::type::value_type) * (bufferCurrIndices[w]-startIdx));
            bufferCnt[w] += (bufferCurrIndices[w]-startIdx);
            uint64_t cnt = bufferCnt[w];
            totalUniqEdgeCnt += cnt;
            std::cerr << "Weight=" << w << " writing " << cnt << " edges.\n";
            outWeightFile[w].seekp(0);
            outWeightFile[w].write(reinterpret_cast<char*>(&cnt), sizeof(uint64_t));
            outWeightFile[w].close();
        }
    }
    std::string sysCommand = "rm -r " + prefix + "tmp*";
    system(sysCommand.c_str());
    logger->info("Duplicated edge count: {}, final unique edge count: {}", tmpedgecnt+colorPairs[0].size(), totalUniqEdgeCnt);
    logger->info("Calculated the weight for all the edges");
    return true;
}

/**
 * Finds Minimum Spanning Forest of color graph using Kruskal Algorithm
 *
 * The algorithm's basic implementation taken from
 * https://www.geeksforgeeks.org/kruskals-minimum-spanning-tree-using-stl-in-c/
 * @return List of connected components in the Minimum Spanning Forest
 */
void MSTMerger::kruskalMSF(AdjList * adjListPtr) {
    uint32_t bucketCnt = numSamples;
    // Create disjoint sets
    DisjointTrees ds(num_colorClasses);
    uint64_t edgeCntr{0}, selectedEdgeCntr{0}, w;

    auto analyze_edge = [&](uint64_t u, uint64_t v, uint64_t w) {
        auto root_of_u = ds.find(u);
        auto root_of_v = ds.find(v);

        // Check if the selected edge is causing a cycle or not
        // (A cycle is induced if u and v belong to the same set)
        if (root_of_u != root_of_v) {
            // Merge two sets
            ds.merge(root_of_u, root_of_v);
            // Current edge will be in the MST
            adjListPtr->storeEdge(u, v, w);
            mstTotalWeight += w;
            selectedEdgeCntr++;
        }
        edgeCntr++;
        if (edgeCntr % 1000000 == 0) {
            std::cerr << "\r" << edgeCntr << " edges processed and "
                      << selectedEdgeCntr << " were selected";
        }
    };
    // Iterate through all sorted edges
    for (uint32_t bucketCntr = 0; bucketCntr < bucketCnt; bucketCntr++) {
        if (bucketCntr < maxWeightInFile) {
            w = bucketCntr + 1;
            TmpFileIterator it(prefix + "w" + std::to_string(w), MAX_ALLOWED_TMP_EDGES_IN_FILE);
            while (not it.end()) {
                auto edge = it.get_val();
                analyze_edge(edge.first, edge.second, w);
                if (selectedEdgeCntr == num_colorClasses) break;
                it.next();
            }
            if (selectedEdgeCntr == num_colorClasses) break;
        } else {
            auto startIdx = ccBitsBucketCnt[bucketCntr];
            auto endIdx = bucketCntr+1==ccBitsBucketCnt.size()?ccBits.size():ccBitsBucketCnt[bucketCntr+1];
            while (startIdx < endIdx) {
                analyze_edge(startIdx, zero, ccBits[startIdx]);
                startIdx++;
            }
        }
    }
    std::cerr << "\n";
    mstTotalWeight++;//1 empty slot for root (zero)
    logger->info("MST Construction finished:"
                 "\n\t# of graph nodes: {}"
                 "\n\t# of graph edges: {}"
                 "\n\t# of merges (mst edges): {}"
                 "\n\tmst weight sum: {}",
                 num_colorClasses, edgeCntr, selectedEdgeCntr, mstTotalWeight);
}

/**
 * calls kruskal algorithm to build an MST of the color graph
 * goes over the MST and fills in the int-vectors parentbv, bbv, and deltabv
 * serializes these three int-vectors as the encoding of color classes
 * @return true if encoding and serializing the DS is successful
 */
bool MSTMerger::encodeColorClassUsingMST() {
    // build mst of color class graph

    logger->info("Running kruskal");
//    usleep(10000000);
    auto adjListPtr = std::make_unique<AdjList>(prefix, num_colorClasses, numSamples);
    kruskalMSF(adjListPtr.get());
    logger->info("After kruskal. Loading the adjacency list.");
//    usleep(10000000);
    adjListPtr->loadCompactedAdjList();
    std::string sysCommand = "rm -r " + prefix + "w*";
    system(sysCommand.c_str());
    sysCommand = "rm -r " + prefix + "tmp*";
    system(sysCommand.c_str());

//    logger->info("after loading the adj list");
//    usleep(10000000);
    uint64_t nodeCntr{0};
    // encode the color classes using mst
    logger->info("Filling ParentBV...");
    sdsl::int_vector<> parentbv(num_colorClasses, 0, ceil(log2(num_colorClasses)));
    sdsl::bit_vector visited(num_colorClasses, 0);
    adjListPtr->hybridTreeWalk(zero, parentbv, visited);
    uint64_t cntr = 0;
    for (auto i = 0; i < visited.size(); i++) {
        if (visited[i] == 0) {
            cntr++;
        }
    }

    std::cerr << "total " << visited.size() << ", visited: " << adjListPtr->visitedCnt << ", not visited: " << cntr << "\n";
    std::vector<uint64_t> thread_deltaOffset_and_parentEnd(nThreads, 0);
    uint64_t idx{0};
    uint64_t bucketSize = std::ceil(parentbv.size() / (double)nThreads);
    uint64_t doubleCheckTotWeight{0};
    for (auto i = 1; i <= adjListPtr->smallerSrcStartIdx.size(); i++) {
        // Limit for adj edges of parent, are the start of the adj edges for next node.
        // So if parent = i-1, limit for adj edges of parent is startIdx[i]
        auto limit = i == adjListPtr->smallerSrcStartIdx.size() ?
                adjListPtr->smallerSrcStartIdx.size() : adjListPtr->smallerSrcStartIdx[i];
        while (idx < limit) {
            auto par = static_cast<uint64_t>(i-1);
            auto val = adjListPtr->smallerSrc[idx];
            uint64_t weight = val & adjListPtr->weightMask;
            uint64_t child = val >> adjListPtr->weightBits;
            if (parentbv[par] == child) {
                std::swap(par, child);
            } else if (parentbv[child] != par) {
                std::cerr << "ERROR! Neither of the two nodes are the parent at index: " << idx << "\n" <<
                             "Expected: " << child << " <-> " << par << "\n"
                             << "Got: " << parentbv[par] << " for " << par <<
                             " and " << parentbv[child] << " for " << child << "\n";
                std::exit(3);
            }
            thread_deltaOffset_and_parentEnd[std::min(static_cast<uint64_t >(nThreads-1), child / bucketSize)] += weight;
            doubleCheckTotWeight += weight;
            idx++;
        }
    }
    if (doubleCheckTotWeight != mstTotalWeight - 1) {
        std::cerr << "ERROR! Weights are not stored properly:\n" <<
        "Expected: " << mstTotalWeight << " Got: " << doubleCheckTotWeight << "\n";
        std::exit(3);
    }
    for (auto i = 1; i < thread_deltaOffset_and_parentEnd.size(); i++) {
        thread_deltaOffset_and_parentEnd[i] += thread_deltaOffset_and_parentEnd[i-1];
    }
    for (auto i = thread_deltaOffset_and_parentEnd.size()-1; i > 0 ; i--) {
        thread_deltaOffset_and_parentEnd[i] = thread_deltaOffset_and_parentEnd[i-1];
//        std::cerr << "thr" << i << " " << thread_deltaOffset_and_parentEnd[i] << "\n";
    }
    thread_deltaOffset_and_parentEnd[0] = 0;

    std::cerr << "\r";
    adjListPtr.reset(nullptr);
    std::cerr << "after deleting MST adjacency list\n";
//    usleep(10000000);

    // fill in deltabv and bbv
    logger->info("Filling DeltaBV and BBV...");
    mst[0].reset(new MSTQuery(prefixes[0], k, k, toBeMergedNumOfSamples[0], logger));
    mst[1].reset(new MSTQuery(prefixes[1], k, k, toBeMergedNumOfSamples[1], logger));
//    mst[0]->loadIdx(prefixes[0]);
//    mst[1]->loadIdx(prefixes[1]);
    sdsl::bit_vector bbv(mstTotalWeight, 0);
    sdsl::int_vector<> deltabv(mstTotalWeight, 0, ceil(log2(numSamples)));
    std::cerr << "after initializing msts in addition to deltabv and bbv\n";
//    usleep(10000000);

//    sdsl::bit_vector::select_1_type sbbv = sdsl::bit_vector::select_1_type(&bbv);
    std::vector<std::thread> threads;
    std::cerr << "Before the start, offsets are:\n";
    for (auto v: thread_deltaOffset_and_parentEnd) {
        std::cerr << v << " ";
    }
    std::cerr << "\n";
    for (uint32_t t = 0; t < nThreads; ++t) {
        threads.emplace_back(std::thread(&MSTMerger::calcDeltasInParallel, this,
                                         t, std::ref(thread_deltaOffset_and_parentEnd[t]),
                                         std::ref(parentbv), std::ref(deltabv), std::ref(bbv), true, adjListPtr.get()));
    }
    for (auto &t : threads) { t.join(); }
    std::cerr << "After calling calcdeltasInParallel, offsets are:\n";
    for (auto v: thread_deltaOffset_and_parentEnd) {
        std::cerr << v << " ";
    }
    std::cerr << "\n";


    bbv[bbv.size()-1] = 1;
    std::cerr << "\r";

    std::cerr << "Done\n";
//    usleep(10000000);

    logger->info("Serializing data structures parentbv, deltabv, & bbv...");
    sdsl::store_to_file(parentbv, std::string(prefix + mantis::PARENTBV_FILE));
    sdsl::store_to_file(deltabv, std::string(prefix + mantis::DELTABV_FILE));
    sdsl::store_to_file(bbv, std::string(prefix + mantis::BOUNDARYBV_FILE));
    logger->info("Done Serializing.");
    return true;
}

void MSTMerger::calcDeltasInParallel(uint32_t threadID, uint64_t &deltaOffset,
                                     sdsl::int_vector<> &parentbv, sdsl::int_vector<> &deltabv,
                                     sdsl::bit_vector &bbv,
                                     bool isMSTBased,
                                     AdjList * adjListPtr) {

    uint64_t deltasKeptInMem = 1000;
    struct Delta {
        uint64_t startingOffset{0};
        std::vector<uint32_t> deltaVals;

        Delta() = default;

        Delta(uint64_t so) {
            startingOffset = so;
        }
    };
    std::vector<Delta> deltas;
    deltas.reserve(deltasKeptInMem);

    colorIdType mst0Zero = ccCnt[0] - 1, mst1Zero = ccCnt[1] - 1;
    auto c1len{static_cast<uint64_t >(ceil(log2(ccCnt[0])))}, c2len{static_cast<uint64_t >(ceil(log2(ccCnt[1])))};
    uint64_t bucketSize = std::ceil(parentbv.size() / (double)nThreads);
    uint64_t s = bucketSize * threadID;;
    uint64_t e = std::min(parentbv.size(), bucketSize * (threadID+1));
    for (uint64_t childIdx = s; childIdx < e; childIdx++) {
//        auto deltaOffset = (childIdx > 0) ? (sbbv(childIdx) + 1) : 0;
        deltas.emplace_back(deltaOffset);
        uint64_t child0{mst0Zero}, parent0{mst0Zero}, child1{mst1Zero}, parent1{mst1Zero};
        if (childIdx != zero) {
            child0 = colorPairs[0][childIdx];
            child1 = colorPairs[1][childIdx];
        }
        if (parentbv[childIdx] != zero) {
            parent0 = colorPairs[0][parentbv[childIdx]];
            parent1 = colorPairs[1][parentbv[childIdx]];
        }

        std::vector<uint32_t> firstDelta, secondDelta;
        if (child0 != parent0)
            firstDelta = getMSTBasedDeltaList(child0, parent0, mst[0].get(), fixed_cache[0],
                                               lru_cache[0][threadID], queryStats[0][threadID]);
        if (child1 != parent1)
            secondDelta = getMSTBasedDeltaList(child1, parent1, mst[1].get(), fixed_cache[1],
                                                lru_cache[1][threadID], queryStats[1][threadID]);
        deltas.back().deltaVals = firstDelta;
        for (auto &v : secondDelta) {
            v += toBeMergedNumOfSamples[0];
            deltas.back().deltaVals.push_back(v);
        }
        deltaOffset += deltas.back().deltaVals.size();
        if (deltas.size() >= deltasKeptInMem) {
            colorMutex.lock();
            for (auto &v : deltas) {
                if (v.startingOffset > 0) {
                    bbv[v.startingOffset - 1] = 1;
                }
                for (auto cntr = 0; cntr < v.deltaVals.size(); cntr++) {
                    deltabv[v.startingOffset + cntr] = v.deltaVals[cntr];
                }
            }
            colorMutex.unlock();
            deltas.clear();
        }
    }

    colorMutex.lock();
    for (auto &v : deltas) {
        if (v.startingOffset > 0) {
            bbv[v.startingOffset - 1] = 1;
        }
        for (auto cntr = 0; cntr < v.deltaVals.size(); cntr++) {
            deltabv[v.startingOffset + cntr] = v.deltaVals[cntr];
        }
    }
    colorMutex.unlock();
}

/**
 * finds the neighbors of each kmer in the cqf,
 * and adds an edge of the element's colorId and its neighbor's (0-based)
 * @param cqf (required to query for existence of neighbors)
 * @param it iterator to the elements of cqf
 */
void MSTMerger::findNeighborEdges(CQF<KeyObject> &cqf, KeyObject &keyobj,
        std::vector<std::pair<uint64_t , uint64_t> > &edgeList,
                                  spp::sparse_hash_map<std::pair<uint64_t , uint64_t >, uint64_t, Custom_Pair_Hasher > & popularEdges) {
    dna::canonical_kmer curr_node(static_cast<int>(k), keyobj.key);
    workItem cur = {curr_node, static_cast<colorIdType>(keyobj.count - 1)};
    uint64_t neighborCnt{0};
    for (auto &nei : neighbors(cqf, cur)) {
        neighborCnt++;
        if (cur.colorId < nei.colorId and popularEdges.find(std::make_pair(cur.colorId, nei.colorId)) == popularEdges.end()) {
//            Edge e(static_cast<colorIdType>(cur.colorId), static_cast<colorIdType>(nei.colorId));

            edgeList.emplace_back(cur.colorId, nei.colorId);
        }
    }
}

/**
 * finds neighbors of a node in cqf
 * @param cqf
 * @param n : work_item containing node and colorId (colorId will be filled)
 * @return set of neighbors for current node n and their colorIds (0-based)
 */
std::set<workItem> MSTMerger::neighbors(CQF<KeyObject> &cqf, workItem n) {
    std::set<workItem> result;
    for (const auto b : dna::bases) {
        uint64_t eqid = 0;
        if (exists(cqf, n.node << b, eqid)) {
            if (eqid != n.colorId)
                result.insert(workItem(n.node << b, eqid));
        }
        if (exists(cqf, b >> n.node, eqid)) {
            if (eqid != n.colorId)
                result.insert(workItem(b >> n.node, eqid));
        }
    }
    return result;
}

/**
 * searches for a kmer in cqf and returns the correct colorId if found
 * which is cqf count value - 1
 * @param cqf
 * @param e : search canonical kmer
 * @param eqid : reference to eqid that'll be set
 * @return true if eqid is found
 */
bool MSTMerger::exists(CQF<KeyObject> &cqf, dna::canonical_kmer e, uint64_t &eqid) {
    KeyObject key(e.val, 0, 0);
    auto eqidtmp = cqf.query(key, QF_NO_LOCK /*QF_KEY_IS_HASH | QF_NO_LOCK*/);
    if (eqidtmp) {
        eqid = eqidtmp - 1;
        return true;
    }
    return false;
}

void MSTMerger::buildMSTBasedColor(uint64_t eqid,
                                   MSTQuery *mst,
                                   LRUCacheMap &lru_cache,
                                   std::vector<uint64_t> &eq,
                                   QueryStats &queryStats,
                                   std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache) {
    RankScores rs(1);

    nonstd::optional<uint64_t> dummy{nonstd::nullopt};

//    auto eq_ptr = lru_cache.lookup_ts(eqid);
//    if (eq_ptr) {
    if ((not fixed_cache.empty()) and fixed_cache.find(eqid) != fixed_cache.end()) {
//        std::cerr << "happens! ";
        eq = fixed_cache[eqid];
        queryStats.cacheCntr++;
//        queryStats.heightDist.push_back(0);
//        queryStats.weightDist.push_back(0);
    } else if (lru_cache.contains(eqid)) {
        eq = lru_cache[eqid];//.get(eqclass_id);
//        eq = *eq_ptr;
        queryStats.cacheCntr++;
//        queryStats.heightDist.push_back(0);
//        queryStats.weightDist.push_back(0);
    } else {
        nonstd::optional<uint64_t> toDecode{nonstd::nullopt};
        queryStats.noCacheCntr++;
        queryStats.trySample = (queryStats.noCacheCntr % 20 == 0);
        toDecode.reset();
        eq = mst->buildColor(eqid, queryStats, &lru_cache, &rs, &fixed_cache, toDecode);
//        auto sp = std::make_shared<std::vector<uint64_t>>(eq);
        lru_cache.emplace(eqid, eq);
//        lru_cache.emplace_ts(eqid, sp);
        if (queryStats.trySample and toDecode and fixed_cache.find(*toDecode) == fixed_cache.end()) {
            auto s = mst->buildColor(*toDecode, queryStats, nullptr, nullptr, &fixed_cache, dummy);
            lru_cache.emplace(*toDecode, s);
//            auto sp1 = std::make_shared<std::vector<uint64_t>>(s);
//            lru_cache.emplace_ts(*toDecode, sp1);
        }
    }
}

uint32_t MSTMerger::mstBasedHammingDist(uint64_t eqid1,
                                        uint64_t eqid2,
                                        MSTQuery *mst,
                                        LRUCacheMap &lru_cache,
                                        QueryStats &queryStats,
                                        std::unordered_map<uint64_t, std::vector<uint64_t>> &fixed_cache) {

    uint32_t dist{0};
    std::vector<uint64_t> eq1, eq2;

    buildMSTBasedColor(eqid1, mst, lru_cache, eq1, queryStats, fixed_cache);
    // fetch the second color ID's BV
    buildMSTBasedColor(eqid2, mst, lru_cache, eq2, queryStats, fixed_cache);
    /// calc distance
    auto i{0}, j{0};
    while (i != eq1.size() and j != eq2.size()) {
        if (eq1[i] == eq2[j]) {
            i++;
            j++;
        } else if (eq1[i] < eq2[j]) {
            i++;
            dist++;
        } else {
            j++;
            dist++;
        }
    }
    if (i != eq1.size()) dist += (eq1.size() - i);
    if (j != eq2.size()) dist += (eq2.size() - j);
    return dist;
}

std::vector<uint32_t> MSTMerger::getMSTBasedDeltaList(uint64_t eqid1, uint64_t eqid2,
                                                      MSTQuery * mstPtr,
                                                      std::unordered_map<uint64_t, std::vector<uint64_t>>& fixed_cache,
                                                      LRUCacheMap &lru_cache,
                                                      QueryStats &queryStats,
                                                      bool verbose) {
    std::vector<uint32_t> res;
    if (eqid1 == eqid2) return res;
    std::vector<uint64_t> eq1, eq2;
    buildMSTBasedColor(eqid1, mstPtr, lru_cache, eq1, queryStats, fixed_cache);
    buildMSTBasedColor(eqid2, mstPtr, lru_cache, eq2, queryStats, fixed_cache);

    if (verbose) {
        std::stringstream ss;
        ss << eqid1 << " child --> ";
        for (auto v: eq1) {
            ss << v << " ";
        }    ss << "\n";
        std::cerr << ss.str();

    }
    if (verbose) {
        std::stringstream ss;
        ss << eqid2 << " parent --> ";
        for (auto v: eq2) {
            ss << v << " ";
        }    ss << "\n";
        std::cerr << ss.str();
    }
    /// calc delta
    auto i{0}, j{0};
    while (i != eq1.size() and j != eq2.size()) {
        if (eq1[i] == eq2[j]) {
            i++;
            j++;
        } else if (eq1[i] < eq2[j]) {
            res.push_back(eq1[i]);
            i++;
        } else {
            res.push_back(eq2[j]);
            j++;
        }
    }
    while (i != eq1.size()) {
        res.push_back(eq1[i]);
        i++;
    }
    while (j != eq2.size()) {
        res.push_back(eq2[j]);
        j++;
    }

    return res; // rely on c++ optimization
}

