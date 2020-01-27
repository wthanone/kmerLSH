#include <chrono>
#include <iostream>
#include <map>
#include <mutex>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <thread>

#include "params.h"
#include "abundance.h"
#include "hash.h"
#include "distance.h"
#include "io.h"

#include "threadpool.hpp"
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

void Cluster(HyperParams& params) {
  vector<Abundance*> unknown_abundance, *unknown_abundance_ptr;
  unknown_abundance_ptr = &unknown_abundance;
  int unknown_abundance_size = 0;
  string head = "";
  int dim = 0;

  auto start_time_total = chrono::high_resolution_clock::now(),
  start_time = chrono::high_resolution_clock::now();

  double max_similarity = 0.9;  // Heuristics for maximum spearman similarity.
  double sim_step = (max_similarity - params.min_similarity)/params.cluster_iteration;
  double threshold = max_similarity;
  int cluster_iteration = 0;
  const int buckets = int(pow(2, params.hash_func_num));

  IO::ReadMatrix(unknown_abundance_ptr, &head, &dim, &unknown_abundance_size, params.scale, params.file_name, params.precision);
  IO::getValue(0);
  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
  cout << "Reading matrix takes secs:\t" << elapsed_read << endl;

  int abundance_per_thread = (unknown_abundance).size() / params.threads_to_use;
  vector<vector<Abundance*>> part_abundance(params.threads_to_use, vector<Abundance*>());
  for (int i = 0; i < params.threads_to_use; ++i) {
	int abundance_to_do = abundance_per_thread;
    if (i == params.threads_to_use - 1) {
      abundance_to_do = (unknown_abundance).size() - i * abundance_per_thread;
	}
	int start_index = i * abundance_per_thread;
	part_abundance[i].insert(part_abundance[i].end(), unknown_abundance.begin()+start_index, unknown_abundance.begin()+start_index+abundance_to_do);
  }

  while (cluster_iteration++ < params.cluster_iteration) {
	vector<vector<Abundance*>> lsh_table(buckets, vector<Abundance*>());
    threshold -= sim_step;
    cout << "Iteration:\t" << cluster_iteration << ", cos sim threshold:\t" << threshold <<" dimension : " <<dim<< endl;

    // Generate LSH.
    hashTable hash_table = LSH::generateHashTable(params.hash_func_num, dim);

    // Apply LSH.
    auto start_time_iteration = chrono::high_resolution_clock::now();
    auto start_time_hashing = chrono::high_resolution_clock::now();

    IO::getValue(1);

	vector<vector<Abundance*>> part_hash_values(params.threads_to_use, vector<Abundance*>());
    vector<vector<int>> part_hash_keys(params.threads_to_use, vector<int>());
	/*
    for (int i = 0; i < params.threads_to_use; ++i) {
      part_hash_values.push_back(new vector<Abundance*>());
      part_hash_keys.push_back(new vector<int>());
    }
	*/
	pool tp(params.threads_to_use);

  	for (int i = 0; i < params.threads_to_use; ++i) {
	  tp.schedule(boost::bind(p_lsh, boost::ref(part_hash_values[i]), boost::ref(part_hash_keys[i]), hash_table, part_abundance[i], buckets));
    }

	tp.wait();
	//tp.clear();
	// Merge thread results.

	IO::getValue(2);
	cout << "Merge thread results " << endl;
	merge_hashtable(&lsh_table, params.hash_func_num, part_hash_values, part_hash_keys);
	IO::getValue(3);
    for (int i = 0; i < params.threads_to_use; ++i) {
      vector<Abundance*>().swap(part_hash_values[i]);
      vector<int>().swap(part_hash_keys[i]);
	  vector<Abundance*>().swap(part_abundance[i]);
    }
    part_hash_keys.clear();
    part_hash_values.clear();

    IO::getValue(2);
    free_abundance(unknown_abundance);
	IO::getValue(3);
    end_time = chrono::high_resolution_clock::now();
    elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time_hashing).count();
    cout << "hashing takes secs:\t" << elapsed_read << endl;

    auto start_time_clustering = chrono::high_resolution_clock::now();
    // Cluster within LSH buckets.

	//vector<vector<Abundance*>> part_abundance(params.threads_to_use, vector<Abundance*>());

	IO::getValue(4);
	int buckets_per_thread = lsh_table.size() / params.threads_to_use;

	//pool tp(params.threads_to_use);

    for (int tid = 0; tid < params.threads_to_use; ++tid) {
      int start_pos = buckets_per_thread * tid;
      int buckets_to_do = buckets_per_thread;

      // If this is the last thread, then take the rest.
      if (tid == params.threads_to_use - 1) {
        buckets_to_do = lsh_table.size() - buckets_per_thread * tid;
      }
	  auto start_each_cluster = chrono::high_resolution_clock::now();
	  cout<<"start p_cluster with " << tid << endl;
      tp.schedule(boost::bind( p_cluster, boost::ref(part_abundance[tid]), &params, &lsh_table, start_pos, start_pos + buckets_to_do, threshold));
	  end_time = chrono::high_resolution_clock::now();
      elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_each_cluster).count();
      cout << "clustering "<< tid << " takes secs:\t" << elapsed_read << endl;
    }

	tp.wait();
	tp.clear();

	for (int i=0; i<buckets; i++){
	  vector<Abundance*>().swap(lsh_table[i]);
	}
	lsh_table.clear();

	//part_abundance.swap(vector<Abundance);
	IO::getValue(6);
	end_time = chrono::high_resolution_clock::now();
	elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time_clustering).count();
	cout << "clustering takes secs:\t" << elapsed_read << endl;
	cout << endl;
  }

  auto start_time_merge = chrono::high_resolution_clock::now();
  merge_abundance(unknown_abundance_ptr, part_abundance);
  IO::getValue(7);
  for(int i=0; i<params.threads_to_use; ++i){
    vector<Abundance*>().swap(part_abundance[i]);
  }
  part_abundance.clear();
  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time_merge).count();
  cout << "merging takes secs:\t" << elapsed_read << endl;

  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time_total).count();
  cout << "msCRUSH algorithm hash+cluster takes (secs): " << elapsed_read << endl;
  cout << "#abundance after clustering: " << (unknown_abundance).size() << endl;

  // Save clusters.
  start_time = chrono::high_resolution_clock::now();
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
  cout << "msCRUSH algorithm in total takes (secs): " << elapsed_read << endl;
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

  auto start_time_total = chrono::high_resolution_clock::now();

  Cluster(params);

  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time_total).count();
  cout << "In all, msCRUSH takes: " << elapsed_read << endl;
  cout << endl;

  return 0;
}
