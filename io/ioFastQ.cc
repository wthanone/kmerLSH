#include "io.h"

namespace Utility {

  void CheckRead(uset_t *g_kmer_ptr, vector<ReadEntry> &read_vec, vector<int> &record_vec, unsigned int num_threads, int tid, float kmer_vote, uset_t *large_kmer_count_ptr) {

  	size_t k = Kmer::k;
  	size_t num_reads = record_vec.size();
  	float kmer_count;
  	float factor;
  	float tmp_ratio;
  	float factor_ratio;
  	double random;
  	uset_t::iterator it;
  	uset_t::iterator it_k;

  	for (size_t i = 0; i < num_reads; i++) {
  		if ( (i % num_threads) == tid ) {
  			ReadEntry read_entry = read_vec[i];
  			if (*(read_entry.s) == '\0') {
  				cout << "\nabnormal read entry skipped\n";
  				cout << read_entry.name << endl;
  				continue;
  			}
  			if (read_entry.len < k+10) continue; //at leaset 10 kmers
  			//cout << read_entry.name << '\t' << read_entry.len << endl;
  			//cout << read_entry.s <<endl;
  			Kmer km(read_entry.s);
  			kmer_count = 0;
  			factor = 0.5;
  			for (size_t j = 0; j <= (read_entry.len - k); j++) {
  				if (j > 0) {
  					km = km.forwardBase((read_entry.s)[j+k-1]);
  				}
  				Kmer tw = km.twin();
  				Kmer rep = (km < tw) ? km : tw;

  				//char tmp[1024];
  				//rep.toString(tmp);
  				//cout << tmp << endl;
  				it = g_kmer_ptr->find(rep);
  				if (it != g_kmer_ptr->end()) {
  					it_k = large_kmer_count_ptr->find(rep);
  					if (it_k != large_kmer_count_ptr->end()){
  						factor = factor + 1;
  						//kmer_count = kmer_count + 0;
  					}
  					else{
  						kmer_count ++;
  					}
  					//char tmp[1024];
  					//km.toString(tmp);
  					//cout << tmp << endl;
  				}
  			}
  			//exit(1);
  			tmp_ratio = kmer_count/(read_entry.len - k);

  			//cout << kmer_count << endl << read_entry.len << endl << k << endl;
  			//cout << tmp_ratio << endl;
  			//exit(1);
  			if (tmp_ratio > kmer_vote) {
  				factor_ratio = factor/kmer_count;
  				random = ((double) rand() / (RAND_MAX));
  				//cout << tmp_ratio << endl;
  				if(random > factor_ratio){
  					record_vec[i] = 1;
  				}
  			}
  		}
  	}
  }
  //extract distinctive reads for each sample
  void ReadExtract(uset_t *g_kmer_ptr, vector<string> &files, string output, float kmer_vote, bool verbose, pool &tp, unsigned int num_threads, uset_t *large_kmer_count_ptr) {

  	//size_t num_reads;
  	FastqFile FQ(files);
  	FILE *of = fopen(output.c_str(), "w");
  	if (of == NULL) {
  		cerr << "Could not open file for writing, " << output << endl;
  		exit(1);
  	}

  	/*pool tp(num_threads);
  	if (verbose) {
  		cout << "start " << num_threads << " threads" << endl;
  	}*/

  	vector<ReadEntry> read_vec(FQ.part_size);
  	size_t part_num_reads;

  	while(FQ.read(read_vec, part_num_reads) > 0) {
  		vector<int> check_vec(part_num_reads, 0);

  		for (unsigned int tid = 0; tid < num_threads; tid++) {
  			tp.schedule(boost::bind(CheckRead, g_kmer_ptr, boost::ref(read_vec), boost::ref(check_vec), num_threads, tid, kmer_vote, large_kmer_count_ptr));
  		}

  		tp.wait(); // wait until all kmers in the read are checked

  		//for (int i = 0; i < part_num_reads; i++) {
  		//	cout << check_vec[i] << "  " << endl;
  		//}
  		//exit(1);
  		int buf_size_factor = 100;
  		char str[8192*buf_size_factor];
  		int tot_write_len = 0;
  		int start_idx = 0;
  		int qual_len;
  		for (size_t i = 0; i < part_num_reads; i++) {
  			if (check_vec[i] == 1) {
  				ReadEntry short_read = read_vec[i];
  				//fprintf(of, "@%s\n", short_read.name);
  				//fprintf(of, "%s\n", short_read.s);
  				//fprintf(of, "+\n");
  				//fprintf(of, "%s\n", short_read.qual);
  				str[start_idx] = '@';
  				memcpy(&str[start_idx + 1], short_read.name, (short_read.name_len)*sizeof(char));
  				str[start_idx + short_read.name_len + 1] = '\n';

  				memcpy(&str[start_idx + short_read.name_len + 2], short_read.s, (short_read.len)*sizeof(char));
  				str[start_idx + short_read.name_len + short_read.len + 2] = '\n';

  				memcpy(&str[start_idx + short_read.name_len + short_read.len + 3], "+\n", 2);

  				qual_len = strlen(short_read.qual);
  				memcpy(&str[start_idx + short_read.name_len + short_read.len + 5], short_read.qual, qual_len*sizeof(char));
  				str[start_idx + short_read.name_len + short_read.len + qual_len + 5] = '\n';

  				tot_write_len += (short_read.name_len + short_read.len + qual_len + 6);

  				start_idx = tot_write_len;
  				if (tot_write_len > 8192*(buf_size_factor-1)) {
  					//cout << "here\n";
  					fwrite(&str, 1, tot_write_len, of);
  					tot_write_len = 0;
  					start_idx = 0;
  				}

  				//fwrite(&str, 1, short_read.name_len + short_read.len + qual_len + 5, of);
  			}
  		}
  		if (tot_write_len > 0 && tot_write_len <= 8192*(buf_size_factor-1))
  			fwrite(&str, 1, tot_write_len, of);
  		check_vec.clear();
  		read_vec.clear();
  		if (part_num_reads < FQ.part_size) break;

  	}
  	fclose(of);
  	FQ.close();
  }

  //added by Mingjie on 09/05/2014
  void GetInput(string input, vector<string> &samples, vector<string> &kmc_names) {
  	ifstream infile (input.c_str());
  	string line;
  	if (infile.is_open()) {
  		while (std::getline(infile, line)) {
  			istringstream ss(line);
  			string sample, kmc_name;
  			ss >> sample >> kmc_name;

  			samples.push_back(sample);
  			kmc_names.push_back(kmc_name);
  			//cout << sample << "\t" << kmc_name << endl;
  		}
  		infile.close();
  	}
  	else cerr << "Unable to open info file";
  }


}  // namespace Utility
