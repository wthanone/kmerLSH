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

namespce Utility{
//identification of distinctive reads using multi-threads
void CheckRead(uset_t *g_kmer_ptr, vector<ReadEntry> &read_vec, vector<int> &record_vec, unsigned int num_threads, int tid, float kmer_vote, uset_t *large_kmer_count_ptr)

//extract distinctive reads for each sample
void ReadExtract(uset_t *g_kmer_ptr, vector<string> &files, string output, float kmer_vote, bool verbose, pool &tp, unsigned int num_threads, uset_t *large_kmer_count_ptr)
}
