#ifndef __CORE_AB_H__
#define __CORE_AB_H__

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

using namespace std;
using namespace Core;

namespace Core {
class AB {
 public:
  static void randAbundance(Abundance* abundance, const Abundance& ab, double scale);
  static void SetAbundance(Abundance* abundance,const vector<string>& kmers, const vector<double>& values, const vector<int>& locs) ;
  static void SetNullAbundance(Abundance* abundance);
  static void convertHTAB(uint16_t **ary_count, vector<uint64_t> &v_kmers, int tot_sample, uint64_t batch_size, streamoff batch_offset, int num_sample, vector<Abundance*>* abVec);
  static void SetConsensus(Abundance* abundance, const Abundance& ab1, const Abundance& ab2);

  static void SetMean(Abundance* abundance,const vector<Abundance*>& local_ab);
  static bool isSameAb(Abundance* ab1, Abundance* ab2);
  static double convert(char const* source, char ** endPtr);
};

}  // namespace Utility
#endif
