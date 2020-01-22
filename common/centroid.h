#ifndef __CORE_CENTROID_H__
#define __CORE_CENTROID_H__

#include <cmath>
#include <iostream>
#include <string>

#include "abundance.h"
#include "params.h"

using std::endl;
using std::fabs;
using std::ostream;
using std::string;

namespace Core {
class Centroid {
 public:
  bool _is_clustered;
  int _num_clustered;
  vector<Abundance*> _genes

  Centroid(): _is_clustered(false), _count(1), _title("default"), _component_titles("default") {};

  ~Abundance() {}

  Abundance(const Abundance& ab) {
    _component_titles = ab._component_titles;
    _count = ab._count;
    _is_clustered = ab._is_clustered;
    _is_consensus = ab._is_consensus;
    _title = ab._title;
    _values = ab._values;
	_locs = ab._locs;
	_num_clustered = ab._num_clustered;
  }

  Abundance& operator=(const Abundance& ab) {
    _component_titles = ab._component_titles;
    _count = ab._count;
    _is_clustered = ab._is_clustered;
    _is_consensus = ab._is_consensus;
    _title = ab._title;
    _values = ab._values;
	_locs = ab._locs;
	_num_clustered = ab._num_clustered;
    return *this;
  }
  
  friend ostream& operator<<(ostream& os, const Abundance& abundance);
/*
  bool shareTopPeaks(const Spectrum& other, double epsilon) {
    int i = 0, j = 0;
    while (i < _top_peak_mz.size() && j < other._top_peak_mz.size()) {
      if (fabs(_top_peak_mz[i] - other._top_peak_mz[j]) < epsilon) {
        return true;
      }
      if (_top_peak_mz[i] < other._top_peak_mz[j]) {
        ++i;
      } else {
        ++j;
      }
    }
    return false;
  }*/
};

inline ostream& operator<<(ostream& os, const Abundance& ab) {
  os << "is_clustered: " << ab._is_clustered << endl;
  os << "is_consensus: " << ab._is_consensus << endl;
  os << "title: " << ab._title << endl;
  os << "component titles: " << ab._component_titles << endl;
  os << "number of clustered: " << ab._num_clustered << endl;
  os << "gene abundances' location and value :" << endl;
  //for (auto value : ab._values) {
  for (int i = 0; i < ab._values.size(); ++i ) {
    os <<ab._locs[i] << ", " << ab._values[i] << endl;
  }
  os << endl;
  return os;
}

}  // namespace Core
#endif
