kmer_result_file = 'query_kmer.res'
lmer_result_file = 'query.res'

def read_result_file(filename, theta):
    result_list = []
    with open(filename, 'r') as f:
        line = f.readline()
        line.strip()
        val = line.split('\t')
        total_kmer = val[1]
        while line:
            line = f.readline()
            if line:
                line.strip()
                val = line.split('\t')
                file_name = val[0]
                count = val[1]
                ratio = int(count)/int(total_kmer)
                if ratio > theta:
                    result_list.append(file_name)
    return result_list

kmer_result = read_result_file(kmer_result_file, 0.0)
print(kmer_result)
lmer_result = read_result_file(lmer_result_file, 0.0)
print(lmer_result)
