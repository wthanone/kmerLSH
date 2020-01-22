#ifndef __CORE_ABUNDANCE_H__
#define __CORE_ABUNDANCE_H__

#include <cmath>
#include <iostream>
#include <string>

#include "params.h"

using std::endl;
using std::fabs;
using std::ostream;
using std::string;

namespace Core {
class Abundance {
 public:
  vector<double> _values;
  vector<int> _locs;
  vector<string> _kmers;

  Abundance() {};

  ~Abundance() {}

  Abundance(const Abundance& ab) {
    _values = ab._values;
	_locs = ab._locs;
	_kmers = ab._kmers;
  }

  Abundance& operator=(const Abundance& ab) {
    _values = ab._values;
	_locs = ab._locs;
	_kmers = ab._kmers;
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
  //os << "number of kmers: " << ab._kmers.size() << endl;
  //os << "gene abundances' location and value :" << endl;
  //for (auto value : ab._values) {
  for (int i = 0; i< ab._kmers.size(); i++){
	os << ab._kmers[i] << "\t" ;
  }
  os << endl;
  for (int i = 0; i < ab._values.size(); ++i ) {
    os <<"(" << ab._locs[i] << ", " << ab._values[i] <<")" ;
  }
  os << endl;
  return os;
}

}  // namespace Core
#endif
