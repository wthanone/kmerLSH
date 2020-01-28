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
class IO {
 public:
  static int parseLine(char* line);
  static int getValue(int i);

  // Use customized char* to double converter, it can support up to precision
  // 0.00001. Changes to [strtof] might be needed to handle higher precision.
  //static void ReadSpectraFromMGF(vector<Spectrum*>* indexed_spectra, 
  //  unordered_map<int, vector<int>>* map_spectra_by_charge,
  //  unordered_map<string, int>* map_ms_titles_to_index,
  //    int* spectra_size, string file_name, double scale, double min_mz,
  //    double max_mz, double precision, int select_topk, int window_mz,
  //    bool peak_normalized, bool remove_precursor, bool verbose=false);
  static void ReadMatrix(vector<Abundance*>* geneAbundances,string* head, int* dim ,int* _size, double scale,  string file_name, double precision);


  // [Deprecated].
  //static void ProcessSpectra(vector<Spectrum>* spectra, int start, int end,
  //    double scale, double min_mz, double precision, int select_topk,
  //    double window_mz);

  // [Deprecated].
  //static void ProcessSpectra(const vector<vector<string>>& chunks, double scale,
  //    double min_mz, double precision, int select_topk, double window_mz);

  //static void SetSpectrum(Spectrum* spectrum,
  //    bool is_clustered, bool is_consensus, int charge,
  //    int count, const EmbededPeaks& embeded_peaks, const Peaks& raw_peaks,
  //    const Peaks& filtered_peaks,
  //    string peptide, double pre_mz, string title, string component_titles,
  //    const vector<double>& top_peak_mz);
  static void randAbundance(Abundance* abundance, const Abundance& ab, double scale);
  static void SetAbundance(Abundance* abundance,const vector<string>& kmers, const vector<double>& values, const vector<int>& locs) ;
  static void SetNullAbundance(Abundance* abundance);
  // Normalize peak intensity with the maximum set to 'scale'.
  //static void Normalize(Peaks* peaks, double scale);

  // Remove adjacent peaks within mz tolerance, keep the strongest peak.
  //static void RemoveAdjacentPeaks(Peaks* peaks, double mz_tolerance);

  // Embed peaks.
  //static void Embed(EmbededPeaks* embeded_peaks, const Peaks& peaks,
   //   double min_mz, double precision, double scale);
  
  // Set consensus spectrum, using just two spectra.
  //static void SetConsensus(Spectrum* consensus, const Spectrum& s1,
  //    const Spectrum& s2, double precision, int topK, double bin_size,
  //    double min_mz, double scale, string title, string component_titles);
  static void SetConsensus(Abundance* abundance, const Abundance& ab1, const Abundance& ab2);

  //static void MergeTwoPeaks(const Peaks& p1, const Peaks& p2, Peaks* peaks);

  // Adjust intensity according to frequency.
  //static void AdaptPeakIntensities(Peaks* peaks, int nSpectra);

  //static void BinTopKPeak(Peaks* top_peaks, const Peaks& peaks, int peaks_size,
  //    int topK, double bin_size);

  // [Deprecated].
  //static vector<double> SelectTopPeakMZ(
  //    const EmbededPeaks& peaks, double precision,int topK = 5);

  //static vector<double> SelectTopPeakMZ(const Peaks& peaks, int topK = 5);
  
  // Customized char* to double converter, only work with precision up to 0.00001
  static void SetMean(Abundance* abundance,const vector<Abundance*>& local_ab);
  static bool isSameAb(Abundance* ab1, Abundance* ab2);
  static double convert(char const* source, char ** endPtr);
};

}  // namespace Utility
#endif
