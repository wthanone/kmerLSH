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
void IO::ReadCluster(vector<Abundance*>* abMatrix, ){


}
void IO::ReadMatrix(vector<Abundance*>* geneAbundances, string* head, int* dim, int* _size, bool normalization,  string file_name ) {
  int sample_cnt = 0;
  int line_cnt = 0;
  double totalgeneAb= 0.0;
  vector<double> nonzero_ab;
  vector<int> nonzero_loc, ids;
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
    nonzero_ab.clear();
    nonzero_loc.clear();
    ids.clear();
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
    ids.push_back((int) geneAb);

    while(lineStart != lineEnd){
      // https://stackoverflow.com/questions/17465061/how-to-parse-space-separated-doubles-in-c-quickly
      lineStart = lineEnd;
      geneAb = strtod(lineStart, &lineEnd);
      ++sample_cnt;
      if (geneAb > 0){
        nonzero_ab.push_back(geneAb);
        if(normalization){totalgeneAb += pow(geneAb, 2);}
        nonzero_loc.push_back(sample_cnt);
      }
    }

    //check using normalization
    if (normalization){
      for(int i = 0; i< nonzero_ab.size(); i++){
        nonzero_ab[i] /= sqrt(totalgeneAb);
      }
    }

    Abundance* abundance = new Abundance();
    SetAbundance(abundance, ids, nonzero_ab, nonzero_loc);
    (*geneAbundances).push_back(abundance);
    ++line_cnt;
    *dim = sample_cnt-1;
    sample_cnt = 0;
    totalgeneAb = 0;
  }
  cout << "read #total lines: " << line_cnt << endl;
  infile.close();
}

void SaveResult(const vector<Abundance*>* abs_all, string out_dir){
  //1) save _ids
  vector<int> ids, locs;
  vector<double> values;
  ofstream out
  out.open(out_dir+"ids.txt");
  for (int i = 0; i < abs_all->size(); ++i) {
		ids = abs_all->at(i)->_ids;
    out << ids[0];
		for(int j=1; j< ids.size(); j++){
			out << "\t" << ids[j] ;
		}
		out<< endl;
	}
  out.close();
  out.open(out_dir+"ids.txt");
  for (int i = 0; i < abs_all->size(); ++i) {
		locs = abs_all->at(i)->_locs;
    out << locs[0];
		for(int j=1; j< locs.size(); j++){
			out << "\t" << locs[j] ;
		}
		out<< endl;
	}
  out.close();
  out.open(out_dir+"ids.txt");
  for (int i = 0; i < abs_all->size(); ++i) {
    values = abs_all->at(i)->_values;
    out << values[0];
    for(int j=1; j< ids.size(); j++){
      out << "\t" << values[j] ;
    }
    out<< endl;
  }
  out.close();
}


void SaveMatrix(const vector<Abundance*>* abs_all, string out_file, string head, int dim) {
	ofstream out(out_file);
	out << head << endl;
	for (int i = 0; i < abs_all->size(); ++i) {
		out << i ;
		int k =0;
		vector<int> locs = abs_all->at(i)->_locs;
		vector<double> values = abs_all->at(i)->_values;
		for(int j=0; j<dim; ++j){
			if (k < locs.size()){
				if(j == locs[k]){
					out << "\t" << values[k];
					k++;
				}
				else{
					out<< "\t0";
				}
			}
			else {
				out<< "\t0";
			}
		}
		out << endl;
	}
	out.close();
}

}  // namespace Utility
