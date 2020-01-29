#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>

#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <time.h>
#include <limits.h>
//#include <atomic>
#include <chrono>

#include <stdint.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <functional>

#include <getopt.h>

#include "Common.hpp"
//#include "utils.hpp"
#include "fastq.hpp"
#include "Kmer.hpp"
#include "HashTables.hpp"


#include "alglib-3.8.2/src/statistics.h"
#include "alglib-3.8.2/src/ap.h"

#include "kmc_api/kmc_file.h"
#include "kmc_reader.hpp"

#include "threadpool.hpp"
#include <boost/bind.hpp>
//#include <boost/thread.hpp>
using namespace std;
//using namespace boost;
using namespace boost::threadpool;
using namespace std::chrono;


//identification of distinctive reads using multi-threads
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


//added by Mingjie on 02/12/2015
//test kmer count in batches
void BatchTest(uint16_t **ary_count, int tot_sample, uint64_t batch_size, streamoff batch_offset, int num_sample1, int num_sample2, vector<uint64_t> &v_kmers, float pvalue_thresh, const vector<Kmer> &kvec, uset_t *g_kmer1_ptr, uset_t *g_kmer2_ptr, uset_t *large_kmer_count_ptr) {
	double freq_array1[num_sample1], freq_array2[num_sample2];
	double bothtails, lefttail, righttail, freq;
	int numZero;
	uint16_t cnt;
	uint64_t tot_cnt;
	uint64_t num_kmers;
	alglib::real_1d_array array_s1, array_s2;

	bothtails = 0;
	lefttail  = 0;
	righttail = 0;

	Kmer km;
	for(uint64_t i=0; i<batch_size; i++) {
		km = kvec[batch_offset+i];
		tot_cnt = 0;
		numZero = 0;
		for(int j=0; j<tot_sample; j++){
			cnt = ary_count[j][i];
			//cout << cnt << endl;
			tot_cnt += cnt;
			if (cnt == 0){
				numZero++;
			}
			num_kmers = v_kmers[j];
			freq = double(cnt)*1000000/num_kmers;

			if(j<num_sample1) {
				freq_array1[j] = freq;
			} else {
				freq_array2[j-num_sample1] = freq;
			}
		}
		//char tmp[1024];
		//km.toString(tmp);
		//cout << tmp << endl;
		//exit(1);
		//if (tot_cnt > tot_sample) {
		array_s1.setcontent(num_sample1, freq_array1);
		array_s2.setcontent(num_sample2, freq_array2);
		alglib::mannwhitneyutest(array_s1, num_sample1, array_s2, num_sample2, bothtails, lefttail, righttail);
		if(numZero < tot_sample * 0.7 ){
			if(lefttail <= pvalue_thresh) {
				g_kmer2_ptr->insert(km);
			} else if (righttail <= pvalue_thresh) {
				g_kmer1_ptr->insert(km);
			}
			if (tot_cnt > tot_sample*2  ){
				large_kmer_count_ptr->insert(km);
			}
		}
	}
}

void buildKHtable(){

    //call kmc
  if (opt.kmc) {
    if (opt.verbose) {
      cout << endl << "...Running KMC using " << opt.num_threads << " threads..." << endl;
      cout << "max_memory: " << opt.max_memory << endl;
      cout << "kmer length: " << opt.k << endl;
    }
    steady_clock::time_point start_time = steady_clock::now();

    for (int i=0; i<tot_sample; i++) {
      string cmd = "kmc -k" + std::to_string(opt.k) + " -r -cs65535 -ci"+std::to_string(opt.count_min) +" -t" + std::to_string(opt.num_threads) + " -m" + std::to_string(opt.max_memory) + " " + samples[i] + " " + kmc_names[i] + " .";
      std::system(cmd.c_str());
    }
    steady_clock::time_point end_time = steady_clock::now();
    auto duration = duration_cast<std::chrono::minutes>(end_time - start_time).count();
    cout << "Running KMC takes " << duration << " minutes\n";
  }
  //exit(0);

  t1 = steady_clock::now();

  if (opt.verbose) {
    cout << endl << "...Reading kmc files..." << endl;
  }

  pool tp(act_threads);

    if (opt.verbose) {
      cout << "start " << tp.size() << " threads" << endl;
    }

    for (int i=0; i<tot_sample; i++) {
      cout << kmc_names[i] << " : start reading" << endl;
      KmcRead(kmc_names[i], kmap_ptr, opt.verbose, tp, act_threads);
      cout << "finishing reading kmc" << endl;
      if (opt.verbose) {
        cout << kmap.size() << endl;
        cout << "load factor: " << kmap.load_factor() << endl;
      }
    }

    //record the kmer in vector based on the order of kmer in the hash table
    //for random access of kmers in BatchTest

    t2 = steady_clock::now();
    auto duration12 = duration_cast<std::chrono::minutes>(t2 - t1).count();
    cout << "Running KmcRead takes " << duration12 << " minutes\n";

    FILE *kmer_file;
    if ((kmer_file = fopen("kmer_set.hex", "w")) == NULL) {
      cerr << "Unable to write kmer_set.hex file\n";
      exit(1);
    }

    kmap_size = kmap.size();
    kvec.resize(kmap_size);

    Kmer km;
    size_t idx_k = 0;
    for(auto it=kmap.cbegin(); !it.is_end();++it) {
      km = it->first;
      kvec[idx_k] = km;
      km.writeBytes(kmer_file);
      idx_k ++;
    }
    fclose(kmer_file);

    t3 = steady_clock::now();
    auto duration23 = duration_cast<std::chrono::minutes>(t3 - t2).count();
    cout << "Creating vector and writing kmer_set.hex takes " << duration23 << " minutes\n";

    //write out kmer count info
    FILE *bin_count_file;
    FILE *log_file; //record total number of kmers in each sample
    if ((bin_count_file = fopen("kmer_count.bin", "wb")) == NULL) {
      cerr << "Unable to write kmer_count.bin file\n";
      exit(1);
    }
    if ((log_file = fopen("kmer_count.log", "w")) == NULL) {
      cerr << "Unable to write kmer_count.log file\n";
      exit(1);
    }

    if (opt.verbose) {
      cout << "\n...Counting kmers and writing to binary..." << endl;
    }

    fprintf(log_file, "%llu", kmap_size);

    for(vector<string>::iterator it=kmc_names.begin(); it!=kmc_names.end(); ++it) {
      string tmp = *it;

      if (opt.verbose) {
      cout << tmp << endl;
      }

      InitializeHT(kmap_ptr);

      kmer_coverage = KmcCount(tmp, kmap_ptr, opt.verbose, tp, act_threads);
      v_kmers.push_back(kmer_coverage);
      fprintf(log_file, "\t%llu", kmer_coverage);

      WriteHT(kmap_ptr, bin_count_file);
    }
    fclose(bin_count_file);
    fclose(log_file);

    t4 = steady_clock::now();
    auto duration34 = duration_cast<std::chrono::minutes>(t4 - t3).count();
    cout << "Running KmcCount and writing to binary takes " << duration34 << " minutes\n";

    //free the memory of kmap
    kmap.clear();
    tp.clear();
}


void SA(int argc, char **argv) {

	SA_ProgramOptions opt;
	SA_ParseOptions(argc,argv,opt);
	ckhmap_t kmap, *kmap_ptr;
	kmap_ptr = &kmap;
  uset_t large_kmer_count, *large_kmer_count_ptr;
  large_kmer_count_ptr = &large_kmer_count;
	size_t kmer_coverage;
	size_t kmap_size;
	std::string filename;
	vector<string> samples1, samples2, samples;
	vector<string> kmc_names1, kmc_names2, kmc_names;
	vector<size_t> v_kmers; //kmer_cov for each sample
	int num_sample1, num_sample2, tot_sample;
	uset_t g_kmer1, g_kmer2, *g_kmer1_ptr, *g_kmer2_ptr;  //differential kmers for two groups in comparison
	g_kmer1_ptr = &g_kmer1;
	g_kmer2_ptr = &g_kmer2;
	vector<Kmer> kvec;

	alglib::real_1d_array array_s1, array_s2;

	int cap_threads = 6;
	steady_clock::time_point t1, t2, t3, t4, t5, t6;

	if (argc < 2) {
		SA_PrintUsage();
		exit(1);
	}

	/*
	if (!SA_CheckOptions(opt)) {
	CountBF_PrintUsage();
	exit(1);
	} */


	// set static global k-value
	Kmer::set_k(opt.k);

	if (opt.verbose) {
		SA_PrintSummary(opt);
	}

	int act_threads;
	act_threads = opt.num_threads > cap_threads ? cap_threads:opt.num_threads;

	GetInput(opt.input1, samples1, kmc_names1);
	GetInput(opt.input2, samples2, kmc_names2);

	samples = samples1;
	samples.insert(samples.end(), samples2.begin(), samples2.end());
	kmc_names = kmc_names1;
	kmc_names.insert(kmc_names.end(), kmc_names2.begin(), kmc_names2.end());

	num_sample1 = samples1.size();
	num_sample2 = samples2.size();
	tot_sample = samples.size();
	if (opt.verbose) {
		cout << endl << "# samples in group 1: " << num_sample1 << endl << "# samples in group 2: " << num_sample2 << endl;
	}
  if (opt.bin){



  }
	else {
		//read in total number of kmers in each sample
		ifstream logStream("kmer_count.log");
		string line;
		getline(logStream, line);
		istringstream ss(line);
		ss >> kmap_size;
		for (int i = 0; i < tot_sample; i++) {
			ss >> kmer_coverage;
			v_kmers.push_back(kmer_coverage);
		}

		//read kmers in kvec to restore the original kmer order in kmap
		ifstream kmer_file("kmer_set.hex");

		kvec.resize(kmap_size);

		uint8_t bytes[(Kmer::MAX_K)/4];
		size_t idx_k = 0;
		for (size_t i = 0; i < kmap_size; i++) {
			kmer_file.read(reinterpret_cast<char*> (&bytes[0]), sizeof(bytes[0])*(Kmer::MAX_K)/4);
			Kmer km(bytes);
			kvec[idx_k] = km;
			idx_k ++;
		}
		t4 = steady_clock::now();
	}

	ifstream inStream("kmer_count.bin", ios::binary);

	const uint64_t batch_thresh = 10000000; //if changed, remember to change 2D array parameter vec_count in readHT and batchTest
	uint64_t batch_size;

	//initialize 2D array
	uint16_t** ary_count = new uint16_t*[tot_sample];
	for (int i=0; i < tot_sample; i++) {
		ary_count[i] = new uint16_t[batch_thresh];
	}

	if (opt.verbose) {
	  cout << "\n...Loading and testing in batches..." << endl;
	}

	//load and test kmer count info in batches
	uint64_t kcnt_rem = kmap_size;
	streamoff batch_offset = 0;
	while(kcnt_rem >= batch_thresh) {
	  batch_size = batch_thresh;

	  ReadHT(inStream, tot_sample, kmap_size, ary_count, batch_size, batch_offset);
	  BatchTest(ary_count, tot_sample, batch_size, batch_offset, num_sample1, num_sample2, v_kmers, opt.pvalue_thresh, kvec, g_kmer1_ptr, g_kmer2_ptr, large_kmer_count_ptr);

	  //cout << g_kmer1.size() << "\t" << g_kmer2.size() <<endl;

	  batch_offset += batch_size;
	  kcnt_rem -= batch_size;
	  if (opt.verbose) {
		  cout << "# loaded kmers: " << batch_offset << endl;
	  }


	}
	if(kcnt_rem > 0 && kcnt_rem < batch_thresh) {
	  batch_size = kcnt_rem;
	  ReadHT(inStream, tot_sample, kmap_size, ary_count, batch_size, batch_offset);
	  BatchTest(ary_count, tot_sample, batch_size, batch_offset, num_sample1, num_sample2, v_kmers, opt.pvalue_thresh, kvec, g_kmer1_ptr, g_kmer2_ptr, large_kmer_count_ptr);
	}

	for (int i = 0; i < tot_sample; ++i) {
		delete [] ary_count[i];
	}
	delete [] ary_count;
	inStream.close();

	t5 = steady_clock::now();
	auto duration45 = duration_cast<std::chrono::minutes>(t5 - t4).count();
	cout << "BatchTesting takes " << duration45 << " minutes\n";

	if (opt.verbose) {
	  cout << "\n# differential kmers in group 1: " << g_kmer1.size() << endl;
	  cout << "# differential kmers in group 2: " << g_kmer2.size() << endl;
	}

	//deallocate the memory of kvec
	vector<Kmer>().swap(kvec);

	if (opt.verbose) {
	  cout << "\n...Writing distinctive reads..."  << endl;
	}


	//unsigned int write_threads;
	//write_threads = opt.num_threads > cap_threads ? opt.num_threads : act_threads;
	pool tp(opt.num_threads);

	if (opt.verbose) {
		cout << "start " << tp.size() << " threads" << endl;
	}

	for(vector<string>::iterator it=samples1.begin(); it!=samples1.end(); ++it) {
	  vector<string> tmp;
	  tmp.push_back(*it);

	  const char * basep = (*it).c_str();
	  std::string base(basename(basep));

	  if (opt.mode == 1) {
		  filename = std::string(opt.output1) + ("_" + base);
		  if (opt.verbose) {
			  cout << "writing to " << filename << endl;
		  }
		  ReadExtract(g_kmer1_ptr, tmp, filename, opt.kmer_vote, opt.verbose, tp, opt.num_threads, large_kmer_count_ptr);

	  } else if (opt.mode == 2) {
		  //down for down-regulated: RNA less than DNA
		  //up for up-regulated: RNA more than DNA
		  filename = "DNA_downreg_" + base;
		  if (opt.verbose) {
			  cout << "writing to " << filename << endl;
		  }
		  ReadExtract(g_kmer1_ptr, tmp, filename, opt.kmer_vote, opt.verbose, tp, opt.num_threads, large_kmer_count_ptr);

		  filename = "DNA_upreg_" + base;
		  if (opt.verbose) {
			  cout << "writing to " << filename << endl;
		  }
		  ReadExtract(g_kmer2_ptr, tmp, filename, opt.kmer_vote, opt.verbose, tp, opt.num_threads, large_kmer_count_ptr);
	  }
	  tmp.clear();
	}

	for(vector<string>::iterator it=samples2.begin(); it!=samples2.end(); ++it) {
	  vector<string> tmp;
	  tmp.push_back(*it);

	  const char * basep = (*it).c_str();
	  std::string base(basename(basep));

	  if (opt.mode == 1) {
		  filename = std::string(opt.output2) + ("_" + base);
		  if (opt.verbose) {
			  cout << "writing to " << filename << endl;
		  }
		  ReadExtract(g_kmer2_ptr, tmp, filename, opt.kmer_vote, opt.verbose, tp, opt.num_threads, large_kmer_count_ptr);

	  } else if (opt.mode == 2) {
		  filename = "RNA_upreg_" + base;
		  if (opt.verbose) {
			  cout << "writing to " << filename << endl;
		  }
		  ReadExtract(g_kmer2_ptr, tmp, filename, opt.kmer_vote, opt.verbose, tp, opt.num_threads, large_kmer_count_ptr);

		  filename = "RNA_downreg_" + base;
		  if (opt.verbose) {
			  cout << "writing to " << filename << endl;
		  }
		  ReadExtract(g_kmer1_ptr, tmp, filename, opt.kmer_vote, opt.verbose, tp, opt.num_threads, large_kmer_count_ptr);
	  }
	  tmp.clear();
	}

	t6 = steady_clock::now();
	auto duration56 = duration_cast<std::chrono::minutes>(t6 - t5).count();
	cout << "ReadExtract takes " << duration56 << " minutes\n";

}
