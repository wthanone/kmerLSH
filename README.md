# kmerLSH


1. Install

> Make

Required packages : KMC3, OpenMP


2. Usage

>./kmerLSH -a file1 -b file2 -o outfile1 -p outfile2
>-a, --input1=STRING             Input filename for metagenome group A
>-b, --input2=STRING             Input filename for metagenome group B

-o, --output1=STRING            Prefix for output of metagenome A

-p, --output2=STRING            Prefix for output of metagenome B

-I, --cluster_iteration=INT           number of iteration for LSH <default 100>

-N, --min_similarity=FLOAT           minimum threshold of similarity <default 0.80>
-K, --kmer_size=INT             Size of k-mers, at most 31
-T, --threads_to_use=INT        Number of threads for running KMC etc. <default 8>
-X, --max-memory=INT            Max memory for running KMC <default 12>
-C, --count-min=INT            Min threshold of k-mer count for running KMC <default 2>
-S, --size_thresh=INT       Threshold of the size of clustering for U-Test <default 500000>
-P, --pval_thresh=FLOAT       For U-test <default 0.01>
-V, --kmer_vote=FLOAT           Percentage threshold of differential k-mers in distinctive reads <default 0.5>
-F, --clust_file_name=STRING           intermediate clustering result file name <default 'clustering_result.txt'>
-M, --mode=STRING                Optional K : run kmc, B : make bin file, C : clustering, E : extract differential reads 
    --verbose                   Print messages during run
    --only                   Run only the setting mode 

