#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>

#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <time.h>
#include <limits.h>
//#include <atomic>
#include <chrono>

#include <stdint.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <functional>

#include <getopt.h>

#include "Common.hpp"
//#include "utils.hpp"
#include "fastq.hpp"
#include "Kmer.hpp"
#include "HashTables.hpp"
#include "abudnace.h"

#include "alglib-3.8.2/src/statistics.h"
#include "alglib-3.8.2/src/ap.h"

#include "kmc_api/kmc_file.h"
#include "kmc_reader.hpp"

#include "threadpool.hpp"
#include <boost/bind.hpp>
//#include <boost/thread.hpp>
using namespace std;
using namespace boost::threadpool;
using namespace std::chrono;

void GetInput(string input, vector<string> &samples, vector<string> &kmc_names);

//initialize hash table with 0
void InitializeHT(ckhmap_t *kmap_ptr);

//added by Mingjie on 02/11/2015
void WriteHT(ckhmap_t *kmap_ptr, FILE *outfile);

//added by Mingjie on 02/12/2015
//load kmer count info in batches
void ReadHT(std::ifstream &infile, int num_sample, uint64_t num_kmer,  uint16_t **ary_count, uint64_t batch_size, streamoff batch_offset);

void buildKHtable(vector<size_t>* v_kmers, pool &tp, bool kmc, bool verbose, size_t ksize, int count_min, unsigned int num_threads, int max_memory, vector<string> samples, vector<string> kmc_names);
