#include "io.h"

namespace Utility {

int IO::parseLine(char* line){
  //This assumes that a digit will be found and the line ends in " Kb".
  int i = strlen(line);
  const char* p = line;
  while (*p <'0' || *p > '9') p++;
  line[i-3] = '\0';
  i = atoi(p);
  return i;
}

int IO::getValue(int i){ //Note: this value is in KB!
  FILE* file = fopen("/proc/self/status", "r");
  int result = -1;
  char line[128];

  while (fgets(line, 128, file) != NULL){
    if (strncmp(line, "VmSize:", 7) == 0){
      result = parseLine(line);
      break;
    }
  }
  fclose(file);
  cout << i << ". Virtual Memory Usage(KB) : " << result << endl;
  return result;
}

double IO::convert(char const* source, char ** endPtr ) {
  char* end;
  int left = strtol( source, &end, 10 );
  double results = left;
  if ( *end == '.' ) {
      char* start = end + 1;
      int right = strtol( start, &end, 10 );
      static double const fracMult[] 
          = { 0.0, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001 };
      results += right * fracMult[ end - start ];
  }
  if ( endPtr != nullptr ) {
      *endPtr = end;
  }
  return results;
}

void IO::ReadMatrix(vector<Abundance*>* geneAbundances, string* head, int* dim, int* _size, double scale,  string file_name, double precisionm) {
  int sample_cnt = 0;
  int gene_cnt=0;
  double geneAb;
  double totalgeneAb= 0.0;
  vector<double> gene_ab;
  vector<int> gene_loc;
  vector<string> kmer_vec;
  string line;
  ifstream infile(file_name);
  if (!infile.is_open()) {
    cout << "Error! file not open!" << endl;
    exit(-1);
  }
   
  cout << "read start : "<< file_name << endl;
  const char* lineStart;
  char* lineEnd;
  
  while (getline(infile, line)){
	gene_ab.clear();
	gene_loc.clear();
	kmer_vec.clear();
	if (line.empty() || line[0] == '#') {
      continue;
    }
	if (line[0] == '\t'){
		*head = line;
		continue;
	}
	//read gene Name
	
	lineStart = line.c_str();
	geneAb = strtod(lineStart, &lineEnd);
	kmer_vec.push_back(to_string((int) geneAb));
	

	while(lineStart != lineEnd){
		// https://stackoverflow.com/questions/17465061/how-to-parse-space-separated-doubles-in-c-quickly
	  
	  lineStart = lineEnd;
	  geneAb = strtod(lineStart, &lineEnd);

	  ++sample_cnt;
	  if (geneAb > 0){
	  	gene_ab.push_back(geneAb * scale);
		totalgeneAb += pow(geneAb * scale, 2);
		gene_loc.push_back(sample_cnt);
	  }  
	}
	//check using normalization
	for(int i = 0; i< gene_ab.size(); i++){
	  gene_ab[i] /= sqrt(totalgeneAb); 
	  //cout << gene_ab[i] << "\t" ;
	}
	//cout << endl;
    Abundance* abundance = new Abundance();
    SetAbundance(abundance, kmer_vec, gene_ab, gene_loc);
		
    (*geneAbundances).push_back(abundance);
    
    ++gene_cnt;
	*dim = sample_cnt-1;
	sample_cnt = 0;
	totalgeneAb = 0;
  }
  cout << "read #total genes: " << gene_cnt << endl;
  
  infile.close();
}

void IO::SetAbundance(Abundance* abundance, const vector<string>& kmers, const vector<double>& values, const vector<int>& locs) {
  (*abundance)._values = values;
  (*abundance)._locs = locs;
  (*abundance)._kmers = kmers;
}
/*
void IO::randAbundance(Abundance* abundance, const Abundance& ab, double scale ){
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
void IO::randomAbundance(Abundance* abundance, int count, double scale){
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
void IO::SetConsensus(Abundance* abundance, const Abundance& ab1, const Abundance& ab2) {
  vector<double> new_ab_value;
  vector<double> ab_value1 = ab1._values, ab_value2 = ab2._values;
  vector<int> new_ab_loc;
  vector<int> ab_loc1 = ab1._locs, ab_loc2 = ab2._locs;
  vector<string> ab1_kmers = ab1._kmers, ab2_kmers = ab2._kmers;
  int cnt = 0, i = 0, j = 0;
  int ab1_count = ab1_kmers.size(), ab2_count = ab2_kmers.size();
  int all_count = ab1_count + ab2_count ;
  ab1_kmers.insert(ab1_kmers.end(),  ab2_kmers.begin(), ab2_kmers.end());
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
  SetAbundance(abundance, ab1_kmers, new_ab_value, new_ab_loc);
}
void IO::SetMean(Abundance* abundance,const vector<Abundance*>& local_ab){
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
bool IO::isSameAb(Abundance* ab1, Abundance* ab2){
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
}  // namespace Utility
