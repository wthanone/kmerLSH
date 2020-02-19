#include "cluster.h"


void p_lsh(vector<Abundance*> &hash_values, vector<int> &hash_keys, hashTable& hash_table, vector<Abundance*> unknown_abs, int buckets) {
  auto start_time = chrono::high_resolution_clock::now();

  vector<Abundance*> local_hash_values;
  vector<int> local_hash_keys;
  // Reserve space to speed.
  const int num_points = unknown_abs.size();
  local_hash_values.reserve(num_points);
  local_hash_keys.reserve(num_points);

  for (int i = 0; i < num_points; ++i) {
    auto& abundance = *(unknown_abs[i]);
    int key = LSH::random_projection(abundance._values,abundance._locs, hash_table);
	//cout << abundance << endl;
	//cout << "key : " << key << endl;
    local_hash_values.push_back(unknown_abs[i]);
    local_hash_keys.push_back(key);
  }
  hash_values.swap(local_hash_values);
  hash_keys.swap(local_hash_keys);
}


void merge_hashtable(vector<vector<Abundance*>>* final_table, int hash_func_num, const vector<vector<Abundance*>>& part_hash_values,  const vector<vector<int>>& part_hash_keys) {
  const int buckets = int(pow(2, hash_func_num));
  vector<vector<Abundance*>> local_table(buckets, vector<Abundance*>());

  for (int i = 0; i < part_hash_keys.size(); ++i) {
    for (int j = 0; j < part_hash_keys[i].size(); ++j) {
      const auto& key = part_hash_keys[i].at(j);
      const auto& ptr = part_hash_values[i].at(j);
      local_table[key].push_back(ptr);
    }
  }
  //cout << "merge hashtable check 1" << endl;
  final_table->swap(local_table);
  //cout << "merge hashtable swap check" << endl;
}

void free_hashtable(vector<vector<Abundance*>>* final_table, int s){
  int size = (final_table->at(s)).size();
  for(int i = 0; i< size; ++i){
	(final_table->at(s)).pop_back();
  }
  vector<Abundance*>().swap(final_table->at(s));
}
void merge_abundance(vector<Abundance*>* final_abundance, const vector<vector<Abundance*>>& part_abundance) {
  vector<Abundance*> local_abundance;
  for (const auto& abundance : part_abundance) {
    local_abundance.insert(local_abundance.end(), abundance.begin(), abundance.end());
  }
  final_abundance->swap(local_abundance);
}
void free_abundance(vector<Abundance*> final_abundance){
  vector<Abundance*> local_abundance;
  for (int i=0; i< final_abundance.size(); ++i){
	(final_abundance).pop_back();
  }
  vector<Abundance*>().swap(final_abundance);
}


void p_cluster(vector<Abundance*> &part_ab, vector<vector<Abundance*>>* lsh_table, int start_pos, int end_pos, double threshold) {
  double distance;
  int size;

  for (int s = start_pos; s < end_pos; ++s) {
    auto& candidates = lsh_table->at(s);
    size = candidates.size();

    auto start_time = chrono::high_resolution_clock::now();
    //cout << "! " << s << "\t" << size << endl;

    for (int i = 1, j = 0 ; i < size;) {
      auto& current = *(candidates[i]);
      //cout << "I " << i << endl;
      for (j = 0; j < i; ++j) {
        auto& candidate = *(candidates[j]);
        distance = Distance::cosine(current._values, current._locs , candidate._values, candidate._locs);
        if (1 - distance >= threshold) {
          Abundance* ab_new = new Abundance();
          AB::SetConsensus(ab_new, current, candidate);
          delete candidates[i];
          delete candidates[j];
          candidates[j] = ab_new;
          candidates[i] = candidates[--size];
          candidates[size] = NULL;
          break;
        }
      }
      if (j == i){
        ++i;
      }
    }
    part_ab.insert(part_ab.end(), candidates.begin(), candidates.begin()+size);
    auto end_time = chrono::high_resolution_clock::now();
    auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
  }
}

void Cluster(string mat_file_name, string result_file_name, double min_similarity, int cluster_iteration, int hash_func_num, unsigned int threads_to_use, bool normalization, bool verbose) {
  vector<Abundance*> unknown_abundance, *unknown_abundance_ptr;
  unknown_abundance_ptr = &unknown_abundance;
  int unknown_abundance_size = 0;
  string head = "";
  int dim = 0;

  auto start_time_total = chrono::high_resolution_clock::now(),
  start_time = chrono::high_resolution_clock::now();

  double max_similarity = 0.9;  // Heuristics for maximum spearman similarity.
  double sim_step = (max_similarity - min_similarity)/cluster_iteration;
  double threshold = max_similarity;
  int iter = 0;
  const int buckets = int(pow(2, hash_func_num));

  IOMat::ReadMatrix(unknown_abundance_ptr, &head, &dim, normalization,  mat_file_name );
  IOMat::ReadCluster(unknown_abundance_ptr, mat_file_name+".clust");

  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
  
  if(verbose){	
	cout << "reading Matrix takes secs:\t" << elapsed_read << endl;
  }
 

  int abundance_per_thread = unknown_abundance.size() / threads_to_use;
  vector<vector<Abundance*>> part_abundance(threads_to_use, vector<Abundance*>());
  for (int i = 0; i < threads_to_use; ++i) {
    int abundance_to_do = abundance_per_thread;
    if (i == threads_to_use - 1) {
      abundance_to_do = unknown_abundance.size() - i * abundance_per_thread;
    }
    int start_index = i * abundance_per_thread;
    part_abundance[i].insert(part_abundance[i].end(), unknown_abundance.begin()+start_index, unknown_abundance.begin()+start_index+abundance_to_do);
  }
  while (iter++ < cluster_iteration) {
    vector<vector<Abundance*>> lsh_table(buckets, vector<Abundance*>());
    threshold -= sim_step;
	if (verbose){
      cout << "Iteration:\t" << iter << ", cos sim threshold:\t" << threshold <<" dimension : " <<dim<< endl;
	}
    // Generate LSH.
    hashTable hash_table = LSH::generateHashTable(hash_func_num, dim);

    // Apply LSH.
    auto start_time_iteration = chrono::high_resolution_clock::now();
    auto start_time_hashing = chrono::high_resolution_clock::now();

    IOMat::getValue(1);

    vector<vector<Abundance*>> part_hash_values(threads_to_use, vector<Abundance*>());
    vector<vector<int>> part_hash_keys(threads_to_use, vector<int>());
    /*
    for (int i = 0; i < params.threads_to_use; ++i) {
      part_hash_values.push_back(new vector<Abundance*>());
      part_hash_keys.push_back(new vector<int>());
    }*/
    pool tp(threads_to_use);

    for (int i = 0; i < threads_to_use; ++i) {
      tp.schedule(boost::bind(p_lsh, boost::ref(part_hash_values[i]), boost::ref(part_hash_keys[i]), hash_table, part_abundance[i], buckets));
    }
    //tp.wait();
    //tp.clear();
    // Merge thread results.

    IOMat::getValue(2);
    cout << "Merge thread results " << endl;

    merge_hashtable(&lsh_table, hash_func_num, part_hash_values, part_hash_keys);

    IOMat::getValue(3);

    for (int i = 0; i < threads_to_use; ++i) {
      vector<Abundance*>().swap(part_hash_values[i]);
      vector<int>().swap(part_hash_keys[i]);
      vector<Abundance*>().swap(part_abundance[i]);
    }
    part_hash_keys.clear();
    part_hash_values.clear();

    IOMat::getValue(2);

    free_abundance(unknown_abundance);

    IOMat::getValue(3);

    end_time = chrono::high_resolution_clock::now();
    elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time_hashing).count();
    cout << "hashing takes secs:\t" << elapsed_read << endl;

    auto start_time_clustering = chrono::high_resolution_clock::now();

    // Cluster within LSH buckets.
    //vector<vector<Abundance*>> part_abundance(params.threads_to_use, vector<Abundance*>());
    IOMat::getValue(4);

    int buckets_per_thread = lsh_table.size() / threads_to_use;

    //pool tp(params.threads_to_use);

    for (int tid = 0; tid < threads_to_use; ++tid) {
      int start_pos = buckets_per_thread * tid;
      int buckets_to_do = buckets_per_thread;

      // If this is the last thread, then take the rest.

      if (tid == threads_to_use - 1) {
        buckets_to_do = lsh_table.size() - buckets_per_thread * tid;
      }

      auto start_each_cluster = chrono::high_resolution_clock::now();
      cout<<"start p_cluster with " << tid << endl;

      tp.schedule(boost::bind( p_cluster, boost::ref(part_abundance[tid]), &lsh_table, start_pos, start_pos + buckets_to_do, threshold));

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
    IOMat::getValue(6);
    end_time = chrono::high_resolution_clock::now();
    elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time_clustering).count();
    cout << "clustering takes secs:\t" << elapsed_read << endl;
    cout << endl;
  }

  auto start_time_merge = chrono::high_resolution_clock::now();
  merge_abundance(unknown_abundance_ptr, part_abundance);

  IOMat::getValue(7);

  for(int i=0; i<threads_to_use; ++i){
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

  IOMat::SaveResult(unknown_abundance_ptr,  result_file_name+".clust", true);
  IOMat::SaveMatrix(unknown_abundance_ptr, result_file_name, dim, true);
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
