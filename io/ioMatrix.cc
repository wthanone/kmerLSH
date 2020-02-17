#include "ioMatrix.h"

namespace Utility {

int IOMat::parseLine(char* line){
  //This assumes that a digit will be found and the line ends in " Kb".
  int i = strlen(line);
  const char* p = line;
  while (*p <'0' || *p > '9') p++;
  line[i-3] = '\0';
  i = atoi(p);
  return i;
}

int IOMat::getValue(int i){ //Note: this value is in KB!
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

double IOMat::convert(char const* source, char ** endPtr ) {
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

void IOMat::ReadCluster(vector<Abundance*>* AbundanceMat, string file_name){
  ifstream infile(file_name+"ids.txt");
  vector<int> new_ids;
  int id;
  int loc = 0;
  string line;

  if (!infile.is_open()) {
    cout << "Error! file not open!" << endl;
    exit(-1);
  }
  cout << "read cluster result : "<< file_name << endl;
  const char* lineStart;
  char* lineEnd;

  while(getline(infile, line)){
    new_ids.clear();

    lineStart = line.c_str();
    id = strtol(lineStart, &lineEnd);
    new_ids.push_back((int) id);

    while(lineStart != lineEnd){
      // https://stackoverflow.com/questions/17465061/how-to-parse-space-separated-doubles-in-c-quickly
      lineStart = lineEnd;
      id = strtol(lineStart, &lineEnd);
      new_ids.push_back((int) id);
    }
    AB::UpdateAbundanceIDs(AbundanceMat->at(loc), new_ids);
    loc++;
  }
}

void IOMat::ReadMatrix(vector<Abundance*>* AbundanceMat, string* head, int* dim, int* _size, bool normalization,  string file_name ) {
  int sample_cnt = 0;
  int line_cnt = 0;
  double totalvalueAb= 0.0;
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
    valueAb = strtod(lineStart, &lineEnd);
    ids.push_back((int) valueAb);

    while(lineStart != lineEnd){
      // https://stackoverflow.com/questions/17465061/how-to-parse-space-separated-doubles-in-c-quickly
      lineStart = lineEnd;
      valueAb = strtod(lineStart, &lineEnd);
      ++sample_cnt;
      if (valueAb > 0){
        nonzero_ab.push_back(valueAb);
        if(normalization){totalvalueAb += pow(valueAb, 2);}
        nonzero_loc.push_back(sample_cnt);
      }
    }

    //check using normalization
    if (normalization){
      for(int i = 0; i< nonzero_ab.size(); i++){
        nonzero_ab[i] /= sqrt(totalvalueAb);
      }
    }

    Abundance* abundance = new Abundance();
    AB::SetAbundance(abundance, ids, nonzero_ab, nonzero_loc);
    (*AbundanceMat).push_back(abundance);
    ++line_cnt;
    *dim = sample_cnt-1;
    sample_cnt = 0;
    totalvalueAb = 0;
  }
  cout << "read #total lines: " << line_cnt << endl;
  infile.close();
}

void IOMat::SaveResult(const vector<Abundance*>* abs_all, string out_file, bool delfile){
  //1) save _ids
  vector<double> values;

  if (delfile){
    if(remove(out_file) != 0){
      perror("File deletion failed");
    }
    else{
      cout << "files are removed" << endl;
    }
  }
  ofstream out
  out.open(out_file, fstream::out | fstream::app);
  for (int i = 0; i < abs_all->size(); ++i) {
		ids = abs_all->at(i)->_ids;
    out << ids[0];
		for(int j=1; j< ids.size(); j++){
			out << "\t" << ids[j] ;
		}
		out<< endl;
	}
  out.close();
}


void IOMat::SaveMatrix(const vector<Abundance*>* abs_all, string out_file, int dim, bool delfile) {
  if (delfile){
    if(remove(out_file) != 0){
      perror("File deletion failed");
    }
    else{
      cout << "files are removed" << endl;
    }
  }
  ofstream out;
  out.open(out_file, fstream::out | fstream::app);

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
