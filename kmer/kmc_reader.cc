
#include "kmc_reader.h"
#include <cstdlib>

void InsertKmer(ckhmap_t *kmap_ptr, vector<Kmer> &kmer_vec, int num_threads, unsigned int tid) {
	size_t num_kmer = kmer_vec.size();
	Kmer km;
	//ckhmap_t::iterator it;
	
	for (size_t i = 0; i < num_kmer; i++) {
		if ( (i % num_threads == tid) ) {
			km = kmer_vec[i];

			Kmer tw = km.twin();
			Kmer rep = (km < tw) ? km : tw;
			
			kmap_ptr->insert(rep, 0);
			//kmap_ptr->insert(KmerIntPair(rep,0));

		}
	}
}
//go through a kmc library, insert kmers to a hash table 
//  and compute total coverages of kmers in the single kmc library
//void KmcRead(string kmc_file_name, ckhmap_t *kmap_ptr, bool verbose, unsigned int num_threads) {
void KmcRead(string kmc_file_name, ckhmap_t *kmap_ptr, bool verbose, unsigned int num_threads) {

	CKMCFile kmer_data_base;
	//int32 i;

	if (!kmer_data_base.OpenForListing(kmc_file_name)) 
		exit(EXIT_FAILURE);
	else {
		uint32 _kmer_length;
		uint32 _mode;
		uint32 _counter_size;
		uint32 _lut_prefix_length;
		uint32 _signature_len;
		uint32 _min_count;
		uint32 _max_count;
		uint64 _total_kmers;

		kmer_data_base.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);

		char str[1024];
		uint32_t cnt;

		vector<Kmer> kmer_vec(_total_kmers);
		
		if (verbose) {
			cout << "total kmers : "<< _total_kmers << endl;
		}

		//TODO: add min_count_to_set and max_count_to_set
		CKmerAPI kmer_object(_kmer_length);

		uint64 idx_k = 0;
		while (kmer_data_base.ReadNextKmer(kmer_object, cnt)) {
			kmer_object.to_string(str);
			//cout << str << endl;
			Kmer km(str);
			
			kmer_vec[idx_k] = km;
			idx_k ++;
		}
		vector<ckhmap_t*> vecKmap_ptr;
		for (int i = 0; i < num_threads; ++i){
			vecKmap_ptr[i].push_back(new ckmap_t ());
		}
		#paragma omp parallel for
		for (unsigned int tid = 0; tid < num_threads; tid++) {
			InsertKmer(vecKmap_ptr[tid], kmer_vec, num_threads, tid);
		}
		
			//Kmer tw = km.twin();
			//Kmer rep = (km < tw) ? km : tw;

			//pair<hmap_t::iterator,bool> ref = kmap_ptr->insert(KmerIntPair(rep,0));
			//kmap_ptr->insert(rep, 0);
			//kmap_ptr->insert(KmerIntPair(rep,0));
			
		kmer_data_base.Close();
		kmer_vec.clear();
	}
}

void SearchKmer(ckhmap_t *kmap_ptr, vector<Kmer> &kmer_vec, vector<uint32> &cnt_vec, int num_threads) {
	size_t num_kmer = kmer_vec.size();
	Kmer km;
	uint32 cnt, tmp;
	
	for (size_t i = 0; i < num_kmer; i++) {
		if ( (i % num_threads == tid) ) {
			km = kmer_vec[i];
			cnt = cnt_vec[i];

			Kmer tw = km.twin();
			Kmer rep = (km < tw) ? km : tw;
			
			if (kmap_ptr->contains(rep)) {
				tmp = kmap_ptr->find(rep) + cnt;
				if (tmp > 65535) tmp = 65535;
				kmap_ptr->update(rep, tmp);
			}
			//it = kmap_ptr->find(rep);
			//if (it != kmap_ptr->end()) {
			//	it->SetVal(it->GetVal()+cnt);
			//}
		}
	}
}

//record kmer counts in hash table for a single sample(kmc library)
size_t KmcCount(string kmc_file_name, ckhmap_t *kmap_ptr, bool verbose, pool &tp, unsigned int num_threads) {

	CKMCFile kmer_data_base;

	if (!kmer_data_base.OpenForListing(kmc_file_name))
		exit(EXIT_FAILURE);
	else {
		uint32 _kmer_length;
		uint32 _mode;
		uint32 _counter_size;
		uint32 _lut_prefix_length;
		uint32 _signature_len;
		uint32 _min_count;
		uint32 _max_count;
		uint64 _total_kmers;

		kmer_data_base.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);

		//float counter;
		char str[1024];
		uint32 cnt;
		size_t tot_coverage = 0;

		vector<Kmer> kmer_vec(_total_kmers);
		vector<uint32> cnt_vec(_total_kmers);
		
		//TODO: add min_count_to_set and max_count_to_set
		CKmerAPI kmer_object(_kmer_length);

		uint64 idx_k = 0;
		while (kmer_data_base.ReadNextKmer(kmer_object, cnt)) {
			kmer_object.to_string(str);
			Kmer km(str);
			
			kmer_vec[idx_k] = km;
			cnt_vec[idx_k] = cnt;
			tot_coverage += cnt;
			idx_k ++;
		}

		/*pool tp(num_threads);
		if (verbose) {
			cout << "start " << num_threads << " threads" << endl;
		}*/
		for (unsigned int tid = 0; tid < num_threads; tid++) {
			SearchKmer( kmap_ptr, boost::ref(kmer_vec), boost::ref(cnt_vec), num_threads, tid);
		}
		tp.wait();
		//tp.clear();

		kmer_data_base.Close();
		kmer_vec.clear();
		cnt_vec.clear();

		return tot_coverage;
	
	}
}
