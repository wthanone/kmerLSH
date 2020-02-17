#ifndef __DISTANCE_H__
#define __DISTANCE_H__

#include <cmath>
#include <utility>
#include <unordered_map>
#include <vector>

#include "../common/params.h"

#include <iostream>
using namespace std;


using std::pair;
using std::unordered_map;
using std::vector;

namespace Core {
class Distance {
 public:
  static double euclidean(const vector<double> lhs, const vector<int> lloc, const vector<double> rhs, const vector<int> rloc);
  static double cosine(const vector<double> lhs, const vector<int> lloc, const vector<double> rhs, const vector<int> rloc);
  //static double spearsman(const GeneAbs& lhs, const GeneAbs& rhs);

};

}  // namespace Core
#endif
