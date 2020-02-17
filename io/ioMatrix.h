#ifndef __UTILITY_IO_MAT__
#define __UTILITY_IO_MAT__

#include <algorithm>
#include <cassert>
#include <cstring>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <map>
#include <random>
#include <queue>
#include <set>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <unistd.h>
#include <utility>
#include <vector>

#include "../common/abundance.h"
#include "../function/funcAB.h"

using namespace std;
using namespace Core;

namespace Utility{
class IOMat {
 public:
   static int parseLine(char* line);

   static int getValue(int i);

   static double convert(char const* source, char ** endPtr ) ;

   static void ReadCluster(vector<Abundance*>* AbundanceMat, string file_name);

   static void ReadMatrix(vector<Abundance*>* AbundanceMat, string* head, int* dim, int* _size, bool normalization,  string file_name ) ;

   static void SaveResult(const vector<Abundance*>* abs_all, string out_dir);

   static void SaveMatrix(const vector<Abundance*>* abs_all, string out_file, string head, int dim);

};

}  // namespace Utility
#endif
