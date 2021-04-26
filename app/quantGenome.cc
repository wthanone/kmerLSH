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
#include "../io/ioFastQ.h"

#include "../utils/fastq.h"

#include "../hash/HashTables.h"
#include "../hash/hash.h"
#include "../hash/lshash.h"

#include "../kmer/kmc_api/kmc_file.h"
#include "../kmer/kmc_reader.h"
#include "../kmer/Kmer.h"

//#include "../utils/threadpool.hpp"
//#include <boost/bind.hpp>

#include "../utils/alglib-3.15.0/src/statistics.h"
#include "../utils/alglib-3.15.0/src/ap.h"
//#include <boost/thread.hpp>

using namespace std;
//using namespace boost;
//using namespace boost::threadpool;
using namespace std::chrono;
using namespace Core;
using namespace Utility;

struct HyperParams
{
  size_t k;
  float min_similarity;
  int samples;
  int max_memory;
  int count_min;
  bool verbose;
  unsigned int threads_to_use;
  string input;
  string output;
  string clust_file_name;
};

void PrintHyperParams(const HyperParams& params) {
  cout << "************ kmers Cluster Params Setting ****************" << endl;
  cout << "group1 result prefix: " << params.output << endl;
  cout << "cluster result file: " << params.clust_file_name << endl;
  cout << "group1 input file: " << params.input << endl;
  cout << "threds to use: " << params.threads_to_use << endl;
  cout << "********************************************************" << endl;
}

void kLSH_PrintUsage() {
  cerr << "kmers LSH "<< KLSH_VERSION << endl << endl;
  cerr << "Clustering of k-mers [from KMC library] " << endl << endl;
  cerr << "Usage: kmerLSH -i1 -i2 -o1 -o2 [options]";
  cerr << endl << endl <<
	"-a, --input1=STRING             Input filename for metagenome group A" << endl <<
	"-o, --output1=STRING            Prefix for output of metagenome A" << endl <<
  "-K, --kmer_size=INT             Size of k-mers, at most " << (int) (Kmer::MAX_K-1)<< endl <<
	"-T, --threads_to_use=INT        Number of threads for running KMC etc. <default 8>" << endl <<
	"-X, --max-memory=INT            Max memory for running KMC <default 12>" << endl <<
	"-C, --count-min=INT            Min threshold of k-mer count for running KMC <default 2>" << endl <<
  "-S, --size_thresh=INT       Threshold of the size of clustering for U-Test <default 500000>" << endl <<
  "-P, --pval_thresh=FLOAT       For U-test <default 0.01>" << endl <<
  "-V, --kmer_vote=FLOAT           Percentage threshold of differential k-mers in distinctive reads <default 0.5>" << endl <<
  "-F, --clust_file_name=STRING           intermediate clustering result file name <default 'clustering_result.txt'>" << endl <<
	"-M, --mode=STRING                Optional K : run kmc, B : make bin file, C : clustering, E : extract differential reads " << endl <<
	"    --verbose                   Print messages during run" << endl << endl <<
	"    --only                   Run only the setting mode " << endl << endl
	;
}

void SetHyperParams(HyperParams* params) {
  (*params).cluster_iteration = 100;  //TODO: Tune this param.
  (*params).min_similarity = 0.80;  //TODO: Tune this param.
  (*params).scale = 1;
  (*params).threads_to_use = 12;
  (*params).max_memory = 12;
  (*params).count_min = 2;
  (*params).k = 23;
  (*params).samples = 10;
  (*params).clust_file_name = "clustering_result.txt";
}

void ParsingCommands(int argc, char*argv[], HyperParams* params) {
  int verbose_flag = 0;
  int only_flag = 0;
  string mode = "";
  const char* opt_string = "o:p:a:b:H:I:N:X:C:T:K:S:P:V:F:M:";
  static struct option long_options[] =
  {
    {"verbose", no_argument,  &verbose_flag, 1},
    {"only", no_argument,  &only_flag, 1},
    {"output", required_argument, 0, 'o'},
  	{"input", required_argument, 0, 'i'},
  	{"max-memory", optional_argument, 0, 'X'},
  	{"count-min", optional_argument, 0, 'C'},
  	{"threads_to_use", optional_argument, 0, 'T'},
    {"kmer_size", optional_argument, 0, 'K'},
    {"clust_file_name", optional_argument, 0, 'F'},
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
          (*params).output = optarg;
          break;
        case 'i':
          (*params).input = optarg;
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
        case 'F':
          (*params).clust_file_name = optarg;
          break;
        default:
          break;
      }
    }
    if (verbose_flag) {
      (*params).verbose = true;
    }
  }

void ReadCluster(vector<Abundance*>* AbundanceMat, int num_samples, string file_name,  streamoff start_line, uint64_t num_lines, bool verbose){
  auto start_time = chrono::high_resolution_clock::now();
  string clust_file_name = file_name +".clust" ;
  uint64_t loc = 0;
  uint64_t id;
  string line;
  const char* lineStart;
  char* lineEnd;
  //vector<int> ids;
  
  /////read binary file
  ifstream inbfile(file_name, ios::in | ios::binary);
  if (!inbfile.is_open()) {
    cout << "Error! file ( " << file_name << " ) not open!" << endl;
    exit(-1);
  }
  inbfile.seekg( start_line * sizeof(float) * num_samples , ios::beg);
  cout << "read start : "<< file_name << ", size :" << num_lines << endl;
  
  vector<vector<float>> values(num_lines, vector<float>(num_samples));
  vector<Abundance*> local_abundance(num_lines);

  for(size_t i =0; i < num_lines; i++){
    //read gene Name
    inbfile.read(reinterpret_cast<char*> (values[i].data()), sizeof(float)*num_samples);
  }
  inbfile.close();

  //////read cluster file
  ifstream infile(clust_file_name);
  if (!infile.is_open()) {
    cout << "Error! file ( " << clust_file_name << " ) not open!" << endl;
    exit(-1);
  }
  cout << "read cluster result : "<< clust_file_name << endl;
  
  for(size_t i =0; i< start_line; i++){
    //getline(infile, line);
	  infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  while(getline(infile, line) && loc < num_lines){

    lineStart = line.c_str();
    id = strtol(lineStart, &lineEnd, 10);
	  //cout << "size : " << id << endl;
	  vector<uint64_t> ids(id);
	  uint64_t numID = 0;

    while(lineStart != lineEnd && numID < ids.size()){
      // https://stackoverflow.com/questions/17465061/how-to-parse-space-separated-floats-in-c-quickly
      //ids.push_back(id);
      lineStart = lineEnd;
      id =strtol(lineStart, &lineEnd, 10);
	  //cout << "numID : " << numID << " id : " << id  <<endl;
	    ids[numID] = id;
	    numID++;
      infile.clear();
    }
    //AB::UpdateAbundanceIDs(AbundanceMat->at(loc), new_ids);
    Abundance* abundance = new Abundance();
    AB::SetAbundance(abundance, ids, values[loc]);
    local_abundance[loc] = abundance;
	  //cout << *abundance << endl;
	  loc++;
	  vector<uint64_t>().swap(ids);
  }
  infile.close();
  swap(*AbundanceMat, local_abundance);
  vector<vector<float>>().swap(values);
  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time).count();
  if(verbose){
    cout << "Reading the result file takes secs: " << elapsed_read << endl;
  }
}

void ReadMatrix(vector<Abundance*>* AbundanceMat, string file_name ) {
  int sample_cnt = 0;
  uint64_t line_cnt = 0;
  float totalvalueAb= 0.0;
  float valueAb;
  vector<uint64_t> ids;
  string line, head;
  ifstream infile(file_name);
  if (!infile.is_open()) {
    cout << "Error! file ( " << file_name << " ) not open!" << endl;
    exit(-1);
  }
  //////read cluster file
  ifstream infile(clust_file_name);
  if (!infile.is_open()) {
    cout << "Error! file ( " << clust_file_name << " ) not open!" << endl;
    exit(-1);
  }
  cout << "read cluster result : "<< clust_file_name << endl;
  
  for(size_t i =0; i< start_line; i++){
    //getline(infile, line);
	  infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  while(getline(infile, line) && loc < num_lines){

    lineStart = line.c_str();
    id = strtol(lineStart, &lineEnd, 10);
	  //cout << "size : " << id << endl;
	  vector<uint64_t> ids(id);
	  uint64_t numID = 0;

    while(lineStart != lineEnd && numID < ids.size()){
      // https://stackoverflow.com/questions/17465061/how-to-parse-space-separated-floats-in-c-quickly
      //ids.push_back(id);
      lineStart = lineEnd;
      id =strtol(lineStart, &lineEnd, 10);
	  //cout << "numID : " << numID << " id : " << id  <<endl;
	    ids[numID] = id;
	    numID++;
      infile.clear();
    }
    //AB::UpdateAbundanceIDs(AbundanceMat->at(loc), new_ids);
    Abundance* abundance = new Abundance();
    AB::SetAbundance(abundance, ids, values[loc]);
    local_abundance[loc] = abundance;
	  //cout << *abundance << endl;
	  loc++;
	  vector<uint64_t>().swap(ids);
  }
  infile.close();
}




void getGenomeQuant(HyperParams& params){
  vector<string> genomes, kmc_names;
  vector<float_t> v_kmers;
  int tot_genome;
  Kmer::set_k(params.k);
  size_t kmap_size;
  float_t kmer_coverage;
  ckhmap_t kmap, *kmap_ptr;
	kmap_ptr = &kmap;

  auto start_time_total = chrono::high_resolution_clock::now();

  GetInput(params.input, genomes, kmc_names);
  
  if(params.verbose){
    cout << "Start to extract the differential reads from raw data" << endl;
  }
  auto start_time_extracting = chrono::high_resolution_clock::now();
  
  string head = "";
  int dim = 0;
  //read clustering result and do statistical testing(WRS)
  
  vector<int> ;

  ifstream logStream("kmer_count.log");
  string line;
  getline(logStream, line);
  istringstream ss(line);
  ss >> kmap_size;
  logStream.close();

  ifstream kmer_file("kmer_set.hex");
	//kvec.resize(kmap_size);
	uint8_t bytes[(Kmer::MAX_K)/4];
	
  for (size_t i = 0; i < kmap_size; i++) {
	  kmer_file.read(reinterpret_cast<char*> (&bytes[0]), sizeof(bytes[0])*(Kmer::MAX_K)/4);
	  Kmer km(bytes);   
	  Kmer tw = km.twin();
		Kmer rep = (km < tw) ? km : tw;	
		kmap.insert(rep, 0);
  }
  readClusterResult(kmap_ptr, params.clust_file_name+".clust");

  for(){
    if (!kmer_data_base.OpenForListing(kmc_file_name)) 
		  exit(EXIT_FAILURE);
	  else {
      vector<uint64_t> 
      CKMCFileInfo info;
		  kmer_data_base.Info(info);
		  char str[ksize];
		  uint32_t cnt;

		  vector<Kmer> kmer_vec(info.total_kmers);
		
		  if (verbose) {
			  cout << "total kmers : "<< info.total_kmers << endl;
		  }

		  //TODO: add min_count_to_set and max_count_to_set
		  CKmerAPI kmer_object(info.kmer_length);
		  uint64 idx_k = 0;
      
		  while (kmer_data_base.ReadNextKmer(kmer_object, cnt)) {
			  kmer_object.to_string(str);
			  //cout << str << endl;
			  Kmer km(str);
			
			  kmer_vec[idx_k] = km;
			  idx_k ++;
		  }

    }
  }
  //clear the id sets and close file
  unordered_set<uint64_t>().swap(g_kmer_id1);
  kmer_file.close();

  // Report time for all the procedures done above.
  auto end_time_extracting = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time_extracting - start_time_extracting).count();
  if(params.verbose){
    cout << "extracting reads takes (secs): " << elapsed_read << endl;
  }

    //clear loaded abundance information
  for(size_t i; i< clusteredab.size(); i++){
    delete clusteredab[i];
  }
  

  // Report time for all the procedures done above.
  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time_total).count();
  cout << "msCRUSH algorithm in total takes (secs): " << elapsed_read << endl;
  
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
