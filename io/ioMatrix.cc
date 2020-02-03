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

void IO::ReadMatrix(vector<Abundance*>* geneAbundances, string* head, int* dim, int* _size, double scale,  string mat_file_name, string tag_file_name, double precisionm) {
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


}  // namespace Utility
