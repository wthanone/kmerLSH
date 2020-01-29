#ifndef __UTILITY_IO_H__
#define __UTILITY_IO_H__

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

#include "abundance.h"
#include "params.h"

using namespace std;
using namespace Core;

namespace Utility {
class AB {
 public:
  static void randAbundance(Abundance* abundance, const Abundance& ab, double scale);
  static void SetAbundance(Abundance* abundance,const vector<string>& kmers, const vector<double>& values, const vector<int>& locs) ;
  static void SetNullAbundance(Abundance* abundance);

  static void SetConsensus(Abundance* abundance, const Abundance& ab1, const Abundance& ab2);

  static void SetMean(Abundance* abundance,const vector<Abundance*>& local_ab);
  static bool isSameAb(Abundance* ab1, Abundance* ab2);
  static double convert(char const* source, char ** endPtr);
};

}  // namespace Utility
#endif
