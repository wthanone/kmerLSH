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
#include <map>
#include <chrono>

#include <stdint.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <functional>

#include <getopt.h>

#include <mutex>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <thread>

#include "../common/Common.h"
#include "../common/abundance.h"

#include "../function/distance.h"
#include "../function/funcAB.h"
#include "../function/cluster.h"

#include "../io/ioMatrix.h"
#include "../io/ioHT.h"

#include "../utils/fastq.h"

#include "../hash/HashTables.h"
#include "../hash/hash.h"
#include "../hash/lshash.h"

#include "../kmer/kmc_api/kmc_file.h"
#include "../kmer/kmc_reader.h"
#include "../kmer/Kmer.h"

#include "../utils/threadpool.hpp"
#include <boost/bind.hpp>

//#include <boost/thread.hpp>

using namespace std;
//using namespace boost;
using namespace boost::threadpool;
using namespace std::chrono;
using namespace Core;
using namespace Utility;

struct HyperParams
{
  size_t k;
  double min_similarity;
  double scale;
  int cluster_iteration;
  int hash_func_num;
  int max_memory;
  int count_min;
  bool verbose;
  bool kmc, bin, mat, clustering;
  unsigned int threads_to_use;
  string kmc_file_name;
  string result_file_name;
  string mat_file_name;
};

void PrintHyperParams(const HyperParams& params) {
  cout << "************ kmers Cluster Params Setting ****************" << endl;
  cout << "cluster iteration: " << params.cluster_iteration << endl;
  cout << "cluster result prefix: " << params.result_file_name << endl;
  cout << "cluster result file path: " << params.mat_file_name << endl;
  cout << "hash function num: " << params.hash_func_num << endl;
  cout << "input file name: " << params.kmc_file_name << endl;
  cout << "min similarity: " << params.min_similarity << endl;
  cout << "threds to use: " << params.threads_to_use << endl;
  cout << "********************************************************" << endl;
}

void kLSH_PrintUsage() {
  cerr << "kmers LSH "<< KLSH_VERSION << endl << endl;
  cerr << "Clustering of k-mers [from KMC library] " << endl << endl;
  cerr << "Usage: kLSH -o -i -m [options]";
  cerr << endl << endl <<
	"-k, --kmer-size=INT             Size of k-mers, at most " << (int) (Kmer::MAX_K-1)<< endl <<
	"-i, --input=STRING             Input filename for metagenome group A" << endl <<
	"-o, --output=STRING            Prefix for output of metagenome A" << endl <<
	"-m, --matrix=STRING            Prefix for output of metagenome A" << endl <<
	"-v, --kmer-vote=FLOAT           Percentage threshold of differential k-mers in distinctive reads <default 0.5>" << endl <<
	"-t, --number-threads=INT        Number of threads for running KMC etc. <default 8>" << endl <<
	"-m, --max-memory=INT            Max memory for running KMC <default 12>" << endl <<
	"-c, --count-min=INT            Min threshold of k-mer count for running KMC <default 2>" << endl <<
	"-M, --mode=STRING                Optional K : run kmc, B : make bin file, M : make matrix file, C : clustering " << endl <<
	"    --verbose                   Print messages during run" << endl << endl <<
	"    --only                   Run only the setting mode " << endl << endl
	;
}

void SetHyperParams(HyperParams* params) {
  (*params).cluster_iteration = 100;  //TODO: Tune this param.
  (*params).hash_func_num = 10;
  (*params).min_similarity = 0.65;  //TODO: Tune this param.
  (*params).scale = 1;
  (*params).threads_to_use = 12;
  (*params).kmc = true;
  (*params).bin = true;
  (*params).mat = true;
  (*params).clustering = true;
  (*params).max_memory = 12;
  (*params).count_min = 2;
  (*params).k = 23;
}

void ParsingCommands(int argc, char*argv[], HyperParams* params) {
  int verbose_flag = 0;
  int only_flag = 0;
  string mode = "";
  const char* opt_string = "o:i:m:H:I:S:X:C:T:K:M:";
  static struct option long_options[] =
  {
    {"verbose", no_argument,  &verbose_flag, 1},
    {"only", no_argument,  &only_flag, 1},
    {"output", required_argument, 0, 'o'},
  	{"input", required_argument, 0, 'i'},
    {"matrix", required_argument, 0, 'm'},
  	{"hash_func_num", optional_argument, 0, 'H'},
    {"cluster_iteration", optional_argument, 0, 'I'},
    {"min_similarity", optional_argument, 0, 'S'},
  	{"max-memory", optional_argument, 0, 'X'},
  	{"count-min", optional_argument, 0, 'C'},
  	{"threads_to_use", optional_argument, 0, 'T'},
    {"kmer-size", optional_argument, 0, 'K'},
    {"mode", optional_argument, 0, 'M'},
  	{0,0,0,0}
  };

    int option_index = 0;
    int c;
    stringstream ss;
    while (true) {
      c = getopt_long(argc,argv,opt_string, long_options, &option_index);

      if (c == -1) {
        break;
      }
      switch (c) {
        case 'o':
          (*params).result_file_name = optarg;
          break;
        case 'i':
          (*params).kmc_file_name = optarg;
  	      break;
        case 'm':
          (*params).mat_file_name = optarg;
  	      break;
        case 'H':
          (*params).hash_func_num = atoi(optarg);
          break;
        case 'I':
          (*params).cluster_iteration = atoi(optarg);
          break;
      	case 'S':
          (*params).min_similarity = atoi(optarg);
          break;
        case 'X':
          (*params).max_memory = atoi(optarg);
          break;
        case 'C':
          (*params).count_min = atoi(optarg);
          break;
        case 'T':
          (*params).threads_to_use = atoi(optarg);
          break;
        case 'K':
          (*params).k = atoi(optarg);
          break;
        case 'M':
          mode = optarg;
          break;
        default:
          break;
      }
    }
    if (verbose_flag) {
      (*params).verbose = true;
    }

    if (only_flag){
      if (mode == "K"){
        (*params).bin = false;
        (*params).mat = false;
        (*params).clustering = false;
      }
      else if (mode == "B"){
        (*params).kmc = false;
        (*params).mat = false;
        (*params).clustering = false;
      }
      else if (mode == "M"){
        (*params).bin = false;
        (*params).kmc = false;
        (*params).clustering = false;
      }
      else if (mode == "C"){
        (*params).bin = false;
        (*params).mat = false;
        (*params).kmc = false;
      }
    } else {
      if (mode == "B"){
        (*params).kmc = false;
      }
      else if (mode == "M"){
        (*params).bin = false;
        (*params).kmc = false;
      }
      else if (mode == "C"){
        (*params).bin = false;
        (*params).mat = false;
        (*params).kmc = false;
      }
    }
  }

void kmerCluster(HyperParams& params){
  vector<string> samples, kmc_names;
  vector<size_t> v_kmers;
  int tot_sample;
  Kmer::set_k(params.k);
  size_t kmap_size, kmer_coverage;

  GetInput(params.kmc_file_name, samples, kmc_names);

  if (params.kmc && params.bin){
	pool tp(params.threads_to_use);
    buildKHtable( &v_kmers,  tp,  params.kmc, params.verbose, params.k, params.count_min, params.threads_to_use, params.max_memory, samples, kmc_names);
  	tp.wait();
	tp.clear();
  }

  if (params.mat){
    //store v_kmer without buildKHtable
    if(!params.bin){
      ifstream logStream("kmer_count.log");
      string line;
      getline(logStream, line);
      istringstream ss(line);
      ss >> kmap_size;
      for (int i = 0; i < tot_sample; i++) {
        ss >> kmer_coverage;
        v_kmers.push_back(kmer_coverage);
      }
    }
    ifstream inStream("kmer_count.bin", ios::binary);

  	const uint64_t batch_thresh = 10000000; //if changed, remember to change 2D array parameter vec_count in readHT and batchTest
  	uint64_t batch_size;

  	//initialize 2D array
  	uint16_t** ary_count = new uint16_t*[tot_sample];
  	for (int i=0; i < tot_sample; i++) {
  		ary_count[i] = new uint16_t[batch_thresh];
  	}

  	if (params.verbose) {
  	  cout << "\n...Loading and testing in batches..." << endl;
  	}
    if(remove(params.mat_file_name.c_str()) != 0){
      perror("File deletion failed");
    }
    else{
      cout << "files are removed" << endl;
    }
  	//load and test kmer count info in batches
  	uint64_t kcnt_rem = kmap_size;
  	streamoff batch_offset = 0;
  	while(kcnt_rem >= batch_thresh) {
  	  batch_size = batch_thresh;

  	  ReadHT(inStream, tot_sample, kmap_size, ary_count, batch_size, batch_offset);
  	  IOMat::convertHTMat(ary_count, v_kmers,  tot_sample, batch_size, batch_offset, params.mat_file_name);

  	  //cout << g_kmer1.size() << "\t" << g_kmer2.size() <<endl;

  	  batch_offset += batch_size;
  	  kcnt_rem -= batch_size;
  	  if (params.verbose) {
  		  cout << "# loaded kmers: " << batch_offset << endl;
  	  }

  	}
  	if(kcnt_rem > 0 && kcnt_rem < batch_thresh) {
  	  batch_size = kcnt_rem;
  	  ReadHT(inStream, tot_sample, kmap_size, ary_count, batch_size, batch_offset);
  	  IOMat::convertHTMat(ary_count, v_kmers,  tot_sample, batch_size, batch_offset, params.mat_file_name);
  	}

  	for (int i = 0; i < tot_sample; ++i) {
  		delete [] ary_count[i];
  	}
  	delete [] ary_count;
  	inStream.close();

  }

  if (params.clustering){
    //read matrix file
    Cluster(params.mat_file_name, params.result_file_name, params.min_similarity, params.cluster_iteration, params.hash_func_num, params.threads_to_use, false, params.verbose);
  }

}

int main(int argc, char **argv) {
	if (argc < 2) {
		kLSH_PrintUsage();
	} else {
    	HyperParams params;
    	SetHyperParams(&params);
    	ParsingCommands(argc, argv, &params);
    	PrintHyperParams(params);
    	kmerCluster(params);
	}
}
