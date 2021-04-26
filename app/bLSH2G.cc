#include <chrono>
#include <iostream>
#include <map>
#include <mutex>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <thread>

#include "../common/params.h"
#include "../common/abundance.h"
#include "../hash/hash.h"
#include "../function/distance.h"
#include "../function/funcAB.h"
#include "../io/ioMatrix.h"

#include "../utils/threadpool.hpp"
#include <boost/bind.hpp>


using namespace boost::threadpool;
using namespace std;
using namespace Core;
using namespace Utility;


struct HyperParams
{
  double min_similarity;
  double scale;
  double precision;
  int cluster_iteration;
  int hash_func_num;
  int rounds_num;
  int threads_to_use;
  string file_name;
  string result_path;
  string result_prefix;
  string cs_path;
};

void PrintHyperParams(const HyperParams& params) {
  cout << "************ msCRUSH Cluster Params Setting ****************" << endl;
  cout << "cluster iteration: " << params.cluster_iteration << endl;
  cout << "cluster result prefix: " << params.result_prefix << endl;
  cout << "cluster result file path: " << params.result_path << endl;
  cout << "hash function num: " << params.hash_func_num << endl;
  cout << "input file name: " << params.file_name << endl;
  cout << "scale : " << params.scale<< endl;
  cout << "min similarity: " << params.min_similarity << endl;
  cout << "precision: " << params.precision << endl;
  cout << "threds to use: " << params.threads_to_use << endl;
  cout << "********************************************************" << endl;
}

void SetHyperParams(HyperParams* params) {
  (*params).cluster_iteration = 100;  //TODO: Tune this param.
  (*params).hash_func_num = 10;
  (*params).min_similarity = 0.65;  //TODO: Tune this param.
  (*params).precision = 0.8;  //TODO: Tune this param.
  (*params).scale = 1;
  (*params).rounds_num = 4;
  (*params).threads_to_use = 20;
}

void SplitCommands(int argc, char*argv[], HyperParams* params) {
  for (int i = 1; i < argc; ++i) {
    if (1 == i) {
      (*params).threads_to_use = stoi(argv[i]);
    } else if (2 == i) {
      (*params).hash_func_num = stoi(argv[i]);
    } else if (3 == i) {
      (*params).cluster_iteration = stoi(argv[i]);
    } else if (4 == i) {
      (*params).min_similarity = stof(argv[i]);
    } else if (5 == i) {
      (*params).result_prefix = string(argv[i]);
	} else if(6 == i){
	  (*params).file_name = string(argv[i]);
    } else {
      cout << "error" << endl;
	  exit(-1);
    }
  }
}

int main (int argc, char *argv[]) {
  if (argc < 7){
    cout << "Missing parameters, at least 9 params." << endl;
    cout << "Usage: ./mscrush_on_general_abundanceMat [threads_to_use] [hash_func_num] [iteration] [min_similarity] [min_mz] [max_mz] [result_prefix] [mgf_file(s)]." << endl;
    return -1;
  }

  HyperParams params;
  SetHyperParams(&params);

  SplitCommands(argc, argv, &params);

  PrintHyperParams(params);

  auto start_time_total = chrono::high_resolution_clock::now(),
  start_time = chrono::high_resolution_clock::now();

  vector<Abundance*> unknown_abundance;

  int unknown_abundance_size = 0;
  string head = "";
  int dim = 0;
  cout << "Loading abundance..." << endl;

  IO::ReadMatrix(&unknown_abundance, &head, &dim, &unknown_abundance_size, params.scale, params.file_name, params.precision);

  IO::getValue(0);

  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
  cout << "Loading spectra takes secs:\t" << elapsed_read << endl;

  cout << "Now we have #unknown profiling: " << unknown_abundance.size() << endl;

  // Run clustering
  Cluster(&unknown_abundance, params, head, dim);

  // Save clusters.
  auto start_time = chrono::high_resolution_clock::now();
  cout << "Saving cluster results starts: " << endl;

  SaveClusters(unknown_abundance_ptr,  params.result_prefix+".out");
  SaveMatrix(unknown_abundance_ptr, params.result_prefix+".mat", head, dim);
  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
  cout << "Save cluster results takes secs: " << elapsed_read << endl;

  // Release allocated memory for spectra w/ charge; save pointers to spectra
  // w/o charge to variable 'spectra_of_no_charge'.
  start_time = chrono::high_resolution_clock::now();

  cout << "Releasing memory starts." << endl;
  for(int i=0; i<(unknown_abundance).size(); ++i){
		delete (unknown_abundance)[i];
  }
  (unknown_abundance).clear();

  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
  cout << "Releasing memory takes: " << elapsed_read << endl;

  // Report time for all the procedures done above.
  elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time_total).count();
  cout << "bLSH algorithm in total takes (secs): " << elapsed_read << endl;

  return 0;
}
