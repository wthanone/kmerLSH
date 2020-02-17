#ifndef __HASH_LSH__
#define __HASH_LSH__

#include <algorithm>
#include <random>
#include <string>
#include <vector>

#include "../common/params.h"

using namespace std;

namespace Hash{
typedef vector<double> hashFunction;
typedef vector<hashFunction> hashTable;
class LSH {
 public:
  static hashFunction generateNormalHashFunc(int size);
  static hashFunction generateUniformHashFunc(int size, double min, double max);
  static hashTable generateHashTable(int htable_size, int hfunction_size);

  // Different LSHs.
  static int random_projection(const vector<double> geneAbs,const vector<int> geneLocs, const hashFunction& function);
  static int random_projection(const vector<double> geneAbs,const vector<int> geneLocs, const hashTable& table);
  static string p_stable(const vector<double> geneAbs, const hashFunction& function, double b, double r);
  static string p_stable(const vector<double> geneAbs, const hashTable& table, double b, double r);
};

}  // namespace Core
#endif
