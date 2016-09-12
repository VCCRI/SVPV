Cluster, merge and genotype samples called with CNVnator to produce VCF file.

**Method**  
1. Calls with sufficiently low e-values are extracted for each sample  
2. All calls are compiled into a unified list  
3. Calls are merged if the have a jaccard index > threshold, with closer calls merged first  
4. Step 3. is repeated until there are no calls left to merge  
5. Samples are genotyped based on the merged calls and the read-depth  
