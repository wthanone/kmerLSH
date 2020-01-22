#include "hash.h"
namespace Core {
  hashFunction LSH::generateNormalHashFunc(int size) {
    hashFunction function(size, 0);

    random_device rd;
    mt19937 gen(rd());

    normal_distribution<> dis(0, 1);
    for (int i = 0; i < size; ++i) {
      function[i] = dis(gen);
    }
    
    //TODO: Might need to consider normalization.
    
    return function;
  }

  //TODO: Test this method.
  hashFunction LSH::generateUniformHashFunc(int size, double min, double max) {
    hashFunction function(size, 0);

    random_device rd;
    mt19937 gen(rd());

    uniform_real_distribution<> dis(min, max);
    for (int i = 0; i < size; ++i) {
      function[i] = dis(gen);
    }
    
    //TODO: Might need to consider normalization.
    
    return function;
  }

  hashTable LSH::generateHashTable(int htable_size, int hfunction_size) {
    hashTable table;
    for (int j = 0; j < htable_size; ++j) {
      table.push_back(generateNormalHashFunc(hfunction_size));
    }
    return table;
  }

  int LSH::random_projection(const vector<double> geneAbs, const vector<int> geneLoc, const hashFunction& function) {
    double sum = 0;
    //for (auto& ab : geneAbs) {
	for (int i = 0; i < geneAbs.size(); ++i){
      sum += function[geneLoc[i]] * geneAbs[i];
    }
    return sum >=0 ? 1 : 0;
  }

  int LSH::random_projection(const vector<double> geneAbs, const vector<int> geneLoc, const hashTable& table) {
    int key = 0;
    for (auto& function : table) {
      key = key * 2 + random_projection(geneAbs, geneLoc, function);
    }
    return key;
  }

  //TODO: Test p_stable functions.
  string LSH::p_stable(const vector<double> geneAbs, const hashFunction& function, double b, double r) {
    double sum = 0;
    for (int i = 0; i < geneAbs.size(); ++i) {
      sum += geneAbs[i] * function[i];
    }
    return to_string(int((sum + b) / r));
  }
  string LSH::p_stable(const vector<double> geneAbs, const hashTable& table, double b, double r) {
    string key;
    for (const auto& function : table) {
      key += p_stable(geneAbs, function, b, r);
    }
    return key;
  }
}

