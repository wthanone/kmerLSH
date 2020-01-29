#include "ht_reader.h"


//initialize hash table with 0
void InitializeHT(ckhmap_t *kmap_ptr) {
	//int i = 0;
	for(auto it=kmap_ptr->begin(); !it.is_end(); ++it) {
		it.set_value(0);
	}
}

//added by Mingjie on 02/11/2015
void WriteHT(ckhmap_t *kmap_ptr, FILE *outfile) {
	//uint32_t kcount;
	uint16_t kcount;
	unsigned int buf_size = 256000;
	uint16_t tmp[buf_size];
	size_t i = 0;
	for(auto it=kmap_ptr->cbegin(); !it.is_end();++it) {
		kcount = it->second;
		//counter_len = Int2PChar((uint64_t)kcount, (uchar*)str);
		//fwrite((const void*) & kcount, sizeof(uint32_t), 1, outfile);
		if (i == buf_size) {
			fwrite(tmp, sizeof(uint16_t), buf_size, outfile);
			i = 0;
		}
		tmp[i] = kcount;
		//fwrite((const void*) & kcount, sizeof(uint16_t), 1, outfile);
		i ++;
	}
	if (i > 0 && i <= buf_size) {
		fwrite(tmp, sizeof(uint16_t), i, outfile);
	}
	
	fflush(outfile);

}

//added by Mingjie on 02/12/2015
//load kmer count info in batches
//revised by Wontack on 08/16/2019
//combine batchTest and ReadHT
void ReadHT(std::ifstream &infile, int num_sample, uint64_t num_kmer, uint64_t batch_size, streamoff batch_offset, vector<uint64_t> &v_kmers, const vector<Kmer> &kvec, vector<Abundance*> *unknwon_abundance) {
	streamoff sample_offset = 0;
	streamoff tot_offset;
	double tot_cnt, freq, num_kmers;	
	vector<double> cntVec, locVec;
	vector<string> kmerVec;
	Km kmer;

	//initialize 2d array
	uint16_t** ary_count = new uint16_t*[num_sample];
	for (int i=0; i < num_sample; i++) {
		ary_count[i] = new uint16_t[batch_thresh];
	}
	
	for(int i=0; i<num_sample; i++) {
		sample_offset = i*num_kmer*sizeof(uint16_t);
		tot_offset = sample_offset + batch_offset*sizeof(uint16_t); 
		if(infile.is_open()) {
			infile.seekg(tot_offset, ios::beg);
			infile.read(reinterpret_cast<char*> (&ary_count[i][0]), sizeof(ary_count[0][0])*batch_size);
		}
	}
	for(uint64_t i=0; i<batch_size; i++){
		tot_cnt = 0;
		cntVec.clear();
		locVec.clear();
		kmerVec.clear();	
		for(int j = 0; j<num_sample; j++){
			num_kmers= double(v_kmers[j]);
			freq = double(ary_count[i][j])/num_kmers;
			if(freq > 0){
				locVec.push_back(j);
				cntVec.push_back(freq);
				tot_cnt += pow(freq, 2);
			}	
		}
		if(tot_cnt > ){
			for(int j=0; j<cntVec.size(); j++){
				cntVec[j] /= sqrt(tot_cnt);
			}
			kmer = kvec[i];
			Abundance* newAb = new Abundance();
			kmerVec.push_back(kmer);
			SetAbundance(abundance, kmerVec, cntVec, locVec);
		}
	}
	
	//exit(1);
}

