Cluster, merge and genotype samples called with CNVnator to produce a VCF file.

**usage**  

Argurment | Description | Notes | Example  
----------|-------------|---------|-------  
-samples | Whitespace delimited file, first column sample names,<br> second column CNVnator calls file | required | 
-o | output VCF file | required | 
-thresh | Jaccard index threshold to use for clustering [default 0.7] | optional | 
-chroms | Comma separated list of chromosomes to process [default all] | optional |  
-e[1-4] | Maximum threshold for given CNVnator e-value [default 1e-6] <br> Note: a call is filtered out in no e-value is below than this value | optional  | -e1 1e-5 -e2 0.003
-header | File containing output vcf file header. <br> If not provided will use supplied hg38 header | optional | 

**Method**  
1. Calls with sufficiently low e-values are extracted for each sample  
2. All calls are compiled into a unified list  
3. Calls are merged if the have a jaccard index > threshold, with closer calls merged first  
4. Step 3. is repeated until there are no calls left to merge  
5. Samples are genotyped based on the merged calls and the read-depth  
