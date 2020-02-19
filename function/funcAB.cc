#include "funcAB.h"

namespace Core {

void AB::SetAbundance(Abundance* abundance, const vector<int>& ids, const vector<double>& values, const vector<int>& locs) {
  (*abundance)._values = values;
  (*abundance)._locs = locs;
  (*abundance)._ids = ids;
}

void AB::UpdateAbundanceIDs(Abundance* abundance, const vector<int>& ids){
  (*abundance)._ids = ids;
}
/*
void AB::randAbundance(Abundance* abundance, const Abundance& ab, double scale ){
  random_device rd;
  mt19937 gen(rd());
  vector<double> new_gene_values;
  vector<double> gene_values = ab._values;
  vector<int> gene_locs = ab._locs;
  string title = ab._titles;
  int size = gene_values.size();
  int count = ab._count;
  double rand = 0.;
  normal_distribution<> dis(0,1);
  for(int i =0; i<size; ++i){
	rand = dis(gen) * scale;
	new_gene_values.push_back(rand + gene_values[i]);

  }
  SetAbundance(abundance, false, count, title, new_gene_values, gene_locs);
}
void AB::randomAbundance(Abundance* abundance, int count, double scale){
  vector<double> gene_ab;
  vector<int> gene_loc;
  double rand;
  string title = "default";
  random_device rd;
  mt19937 gen(rd());

  normal_distribution<> dis(0,1);
  for (int i=0; i<count; ++i){
	rand = dis(gen) * scale;
	if (rand > 0){
	  gene_ab.push_back(rand);
	  gene_loc.push_back(i);
	}
  }
  SetAbundance(abundance, false, false, count,1, title, title, gene_ab, gene_loc);
}
*/
void AB::SetConsensus(Abundance* abundance, const Abundance& ab1, const Abundance& ab2) {
  vector<double> new_ab_value;
  vector<double> ab_value1 = ab1._values, ab_value2 = ab2._values;
  vector<int> new_ab_loc;
  vector<int> ab_loc1 = ab1._locs, ab_loc2 = ab2._locs;
  vector<int> ab1_ids = ab1._ids, ab2_ids = ab2._ids;
  int cnt = 0, i = 0, j = 0;
  int ab1_count = ab1_ids.size(), ab2_count = ab2_ids.size();
  int all_count = ab1_count + ab2_count ;
  ab1_ids.insert(ab1_ids.end(),  ab2_ids.begin(), ab2_ids.end());

  while(i < ab_loc1.size() && j < ab_loc2.size()){
    if (ab_loc1[i] == ab_loc2[j]){
      new_ab_value.push_back((ab_value1[i]*ab1_count+ab_value2[i]*ab2_count)/all_count);
      new_ab_loc.push_back(ab_loc1[i]);
      ++j;
      ++i;
    }
    else if(ab_loc1[i] < ab_loc2[j]) {
      new_ab_value.push_back(ab_value1[i]*ab1_count/all_count);
      new_ab_loc.push_back(ab_loc1[i]);
      ++i;
    }
    else {
      new_ab_value.push_back(ab_value2[j]*ab2_count/all_count);
      new_ab_loc.push_back(ab_loc2[j]);
      ++j;
    }
    ++cnt;
  }

  while (i < ab_loc1.size()) {
    new_ab_value.push_back(ab_value1[i]*ab1_count/all_count);
    new_ab_loc.push_back(ab_loc1[i]);
    ++i;
    ++cnt;
  }
  while (j < ab_loc2.size()) {
    new_ab_value.push_back(ab_value2[j]*ab2_count/all_count);
    new_ab_loc.push_back(ab_loc2[j]);
    ++j;
    ++cnt;
  }

  SetAbundance(abundance, ab1_ids, new_ab_value, new_ab_loc);
}

void AB::SetMean(Abundance* abundance,const vector<Abundance*>& local_ab){
  int length = local_ab.size();
  Abundance* tmp = new Abundance();
  auto& previous = *(local_ab[0]);
  for(int i =1; i < length; ++i){
    //cout << "SetMean : " << i << endl;
    auto& current = *(local_ab[i]);
    SetConsensus(tmp, current, previous);
    previous = *tmp;
  }
  abundance = tmp;
}

bool AB::isSameAb(Abundance* ab1, Abundance* ab2){
  bool out = false;
  vector<double> ab_values1 = ab1->_values;
  vector<double> ab_values2 = ab2->_values;
  vector<int> ab_locs1 = ab1->_locs;
  vector<int> ab_locs2 = ab2->_locs;
  cout << "isSameAb start" << endl;
  if(ab_locs1.size() != ab_locs2.size()){
    cout << "isSame not same size " << ab_locs1.size() << " " << ab_locs2.size() << endl;
    return out;
  }
  else {
    cout << "isSameAb same size" << endl;
    for(int i=0; i<ab_locs1.size(); ++i){
      if (ab_locs1[i] != ab_locs2[i] || ab_values1[i] != ab_values2[i]){
        return out;
      }
    }
  }
  out = true;
  return out;
}
}
