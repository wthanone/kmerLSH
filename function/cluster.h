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
#include "../hash/lshash.h"
#include "distance.h"
#include "funcAB.h"

#include "threadpool.hpp"
#include <boost/bind.hpp>


using namespace boost::threadpool;
using namespace std;
using namespace Core;
using namespace Utility;

void p_lsh(vector<Abundance*> &hash_values, vector<int> &hash_keys, hashTable& hash_table, vector<Abundance*> unknown_abs, int buckets);

void merge_hashtable(vector<vector<Abundance*>>* final_table, int hash_func_num, const vector<vector<Abundance*>>& part_hash_values,  const vector<vector<int>>& part_hash_keys);

void free_hashtable(vector<vector<Abundance*>>* final_table, int s);


void p_cluster(vector<Abundance*> &part_ab, vector<vector<Abundance*>>* lsh_table, int start_pos, int end_pos, double threshold) ;

void Cluster(vector<Abundance*>* unknown_abundance, double min_similarity, int cluster_iteration, int hash_func_num, unsigned int threads_to_use, bool verbose) ;
