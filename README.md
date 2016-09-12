Structural Variant Prediction Viewer  
------------------------------------

###Requirements
* Python 2.7
* R
* [SAMtools and BCFtools](https://github.com/samtools) (version 1.3)
* Linux environment, or access to linux via ssh (X11 required for GUI)
* ImageMagick (GUI only, for viewing pdfs)

###Installation
* clone this repository  
`git clone https://github.com/VCCRI/SVPV.git SVPV`
* test everything is working  
`python /SVPV/SVPV.py -example`
* test the gui is working  
`python /SVPV/SVPV.py -example -gui`
* all done!  

**Note:** SAMtools and BCFtools must be executable by typing 'samtools' and 'bcftools' into the terminal.

###Usage

**example:**  
```
python SVPV.py -vcf input_svs.vcf -samples sample1,sample2 -aln alignment1.bam,alignment2.sam -o /out/directory/
```

Run args: | Description | Notes
----------|-------------|------
-vcf | Primary structural variant prediction vcf/bcf file | required
-o | Output directory | required
-aln | Comma separated list of alignment files | required*
-samples | Comma separated list of samples to view, names must be the same as in vcf | required*
-gui | run in gui mode. | optional
-alt_vcf | Alternate structural variant prediction vcf/bcf file, <br> called on the same set of samples as primary | optional
-ref_vcf | Reference structural variant vcf/bcf file for annotation | optional
-ref_gene | Refseq genes file for annotation | optional
-manifest | Whitespace delimited file, first column sample names, <br> second column alignment file path. Overrides '-samples' and '-aln' if also given. | optional
*'-samples' and '-aln' not required if '-manifest' is supplied.  

Filter args: | Description | Example  
-------------|-------------|---------
-max_len | maximum length of structural variants (bp)
-min_len | minimum length of structural variants (bp)
-af | Allele frequency threshold |  -af <0.1  <br> for SV with frequency less than 0.1.
-gts | Specify genotypes of given samples | sample1:0/1,1/1;sample2:1/1
-chrom | Restrict to comma separated list of chromosomes
-svtype | Restrict to given SV type (DEL/DUP/CNV/INV)
-rgi | Restrict to SVs that intersect refGenes, <br>'-ref_gene' must be supplied
-exonic | Restrict to SVs that intersect exons of refGenes, <br>'-ref_gene' must be supplied

Plot args: | Default | Description
-----------|---------|------------
-d | 1 | force sequencing depth plot on or off
-or | 1 | force orphaned reads plot on or off
-v | 1 | force inverted pairs plot on or off
-ss | 1 | force same strand pairs plot on or off
-hc | 1 | force hardclipped reads plot on or off
-se | 0 | force SAM 'secondary alignment' plot on or off
-su | 0 | force SAM 'supplementary alignment' plot on or off
-i | 1 | force inferred insert size plot on or off
-r | 1 | force refgenes plot on or off
-af | 1 | force allele frequency plot on or off
-l | 1 | force plot legend on or off

**extended example:**  
```
python SVPV.py -vcf caller1_svs.vcf -samples sample1,sample2,sample3 -aln alignment1.bam,alignment2.bam,alignment3.bam -o /out/directory/ -alt_vcf caller2_svs.vcf -ref_vcf 1000_genomes_svs.vcf -ref_gene hg38.refgene.txt -max_len 100000 -af <0.25 -gts sample1:1/1,0/1;sample3:0/0 -svtype DEL -exonic -ss 0 -se 1
```
