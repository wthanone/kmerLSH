#include "distance.h"
namespace Core{
/*
double Distance::spearsman(const Abundance& lhs, const Abundance& rhs){
   double distance = 0.
   alglib::real_1d_array array_s1, array_s2;
   arr
   array_s1.setcontent(num_sample1, freq_array1);
   array_s2.setcontent(num_sample2, freq_array2);
   alglib::mannwhitneyutest(array_s1, num_sample1, array_s2, num_sample2, bothtails, lefttail, righttail);

   //distance = 


}
*/
//TODO: Implement a highly efficient method like 'cosine'.

double Distance::euclidean(const vector<double> lhs,const vector<int> lloc,  const vector<double> rhs, const vector<int> rloc) {
  double distance = 0;
  int i = 0, j = 0;
  while (i < lhs.size() && j < rhs.size() ) {
    if (lloc[i] == rloc[j]){
      distance += pow( rhs[i] - lhs[i] ,2);
      ++i;
      ++j;
    }
    else if (lloc[i] < rloc[j]) {
      distance += lhs[i] * lhs[i];
      ++i;
    }
    else {
      distance += rhs[j] * rhs[j] ;
      ++j;
    }
  }
  while (i < lhs.size()){
    distance += lhs[i] * lhs[i];
    ++i;
  }
  while (j < rhs.size()){
    distance += rhs[j] * rhs[j];
    ++j;
  }
  return sqrt(distance);
}

double Distance::cosine(const vector<double> lhs,const vector<int> lloc,  const vector<double> rhs, const vector<int> rloc) {
  double similarity = 0;
  double magnitude_lhs = 0, magnitude_rhs = 0;
  int i = 0, j = 0;
  //cout << "cosine lhs : "<< lhs.size() << " rhs : " << rhs.size() << endl;
  while (i < lhs.size() && j < rhs.size() ) {
	if (lloc[i] == rloc[j]){
      similarity += lhs[i] * rhs[i];
      magnitude_lhs += lhs[i] * lhs[i]; 
      magnitude_rhs += rhs[i] * rhs[i];
	  ++i;
	  ++j;
	}
	else if (lloc[i] < rloc[j]) {
	  magnitude_lhs += lhs[i] * lhs[i];
	  ++i;
	}
	else {
	  magnitude_rhs += rhs[j] * rhs[j] ;
	  ++j;
	}
  }
  while (i < lhs.size()){
	magnitude_lhs += lhs[i] * lhs[i];
	++i;
  }
  while (j < rhs.size()){
	magnitude_rhs += rhs[j] * rhs[j];
	++j;
  }
  similarity /= sqrt(magnitude_lhs) * sqrt(magnitude_rhs);
  return 1 - similarity;
}
  


}  //namespace Core
