#include <fstream>
#include "kmer.h"
#include "minimizer.h"

/*return the integer representation of the base */
inline char Kmer::map_int(uint8_t base)
{
	switch(base) {
		case DNA_MAP::A: { return 'A'; }
		case DNA_MAP::T: { return 'T'; }
		case DNA_MAP::C: { return 'C'; }
		case DNA_MAP::G: { return 'G'; }
		default:  { return DNA_MAP::G+1; }
	}
}



/**
 * Converts a string of "ATCG" to a uint64_t
 * where each character is represented by using only two bits
 */
uint64_t Kmer::str_to_int(string str)
{
	uint64_t strint = 0;
	for (auto it = str.begin(); it != str.end(); it++) {
		uint8_t curr = 0;
		switch (*it) {
			case 'A': { curr = DNA_MAP::A; break; }
			case 'T': { curr = DNA_MAP::T; break; }
			case 'C': { curr = DNA_MAP::C; break; }
			case 'G': { curr = DNA_MAP::G; break; }
		}
		strint = strint | curr;
		strint = strint << 2;
	}
	return strint >> 2;
}

std::string Kmer::int_to_str_lmer(__int128_t kmer, uint64_t kmer_size)
{

    char BASES[] = {'A', 'C', 'G', 'T'};
    uint8_t base;
    std::string str;
    for (uint32_t i = kmer_size; i > 0; i--) {
        base = (kmer >> (i*2-2)) & 3ULL;
        char chr = BASES[base];
        str.push_back(chr);
    }
    return str;
}

/**
 * Converts a uint64_t to a string of "ACTG"
 * where each character is represented by using only two bits
 */
string Kmer::int_to_str(uint64_t kmer, uint64_t kmer_size)
{
	uint8_t base;
	string str;
	for (int i=kmer_size; i>0; i--) {
		base = (kmer >> (i*2-2)) & 3ULL;
		char chr = Kmer::map_int(base);
		str.push_back(chr);
	}
	return str;
}


/* Calculate the revsese complement of a kmer */
__int128_t Kmer::reverse_complement(__int128_t kmer, uint64_t kmer_size)
{
	__int128_t rc = 0;
	uint8_t base = 0;
	for (uint32_t i = 0; i < kmer_size; i++) {
		base = kmer & 3ULL;
		base = reverse_complement_base(base);
		kmer >>= 2;
		rc |= base;
		rc <<= 2;
	}
	rc >>=2;
	return rc;
}

/* Compare the kmer and its reverse complement and return the result 
 * Return true if the kmer is greater than or equal to its
 * reverse complement. 
 * */
bool Kmer::compare_kmers(__int128_t kmer, __int128_t kmer_rev)
{
	return kmer >= kmer_rev;
}

mantis::QuerySets Kmer::parse_kmers(const char *filename, uint64_t kmer_size,
																		uint64_t& total_kmers,
																		bool is_bulk,
									//nonstd::optional<std::unordered_map<mantis::KmerHash, uint64_t>> &uniqueKmers
									std::unordered_map<mantis::KmerHash, uint64_t> &uniqueKmers, uint64_t lmer_size) {
    MinimizerScanner scanner(kmer_size, lmer_size, 0, true, 0);
    mantis::QuerySets multi_kmers;
	total_kmers = 0;
	std::ifstream ipfile(filename);
	std::string read;
	while (ipfile >> read) {
		mantis::QuerySet kmers_set;

start_read:
		if (read.length() < kmer_size)
			continue;
		{
            scanner.LoadSequence(read);
            uint64_t *mmp;
            while ((mmp = scanner.NextMinimizer()) != nullptr){
                kmers_set.insert(*mmp);
                if (is_bulk)
                    if (uniqueKmers.find(*mmp) == uniqueKmers.end())
                        uniqueKmers[*mmp] = 0;
            }
		}
next_read:
		total_kmers += kmers_set.size();
		//if (kmers_set.size() != kmers.size())
		//std::cout << "set size: " << kmers_set.size() << " vector size: " << kmers.size() << endl;
		multi_kmers.push_back(kmers_set);
	}
	return multi_kmers;
}

// A=1, C=0, T=2, G=3
std::string Kmer::generate_random_string(uint64_t len) {
	std::string transcript;
	for (uint32_t i = 0; i < len; i++) {
		uint8_t c = rand() % 4;
		transcript += Kmer::map_int(c);
	}
	return transcript;
}
