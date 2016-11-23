Structural Variant Prediction Viewer  
------------------------------------
View predicted structural variant regions from paired-end whole genome sequencing alignments and compare calls from differenct structural variant prediction algorithms. Statistics related to structural variants are presented in a form that allows users to visually identify false postive calls. Input is a set alignment files (SAM/BAM/CRAM format), and a set of structural variant predictions on these alignments (VCF files). Output is a set of pdf files of structural variant regions.

Deletions (DEL), duplications (DUP), copy number variations (CNV) and inversions (INV) are supported. Translocations (TRA) and insertions (INS) are not supported.

###Requirements
**Basic**  
* Python 2.7 and NumPy
* R
* [SAMtools and BCFtools](https://github.com/samtools) (version 1.3)
* Linux environment, or access to linux via ssh  
  
**GUI**  
* X11 if running over ssh
* python 2.7 tkinter
* [GraphicsMagick](http://www.graphicsmagick.org/) or ImageMagick ('display')
  * Recomended: GraphicsMagick build w/ --enable-magick-compat
  * or other X11 pdf viewer specified in 'config.py'

###Installation
* navigate to desired install directory and clone this repository  
`git clone https://github.com/VCCRI/SVPV.git SVPV`
* test everything is working  
`python ./SVPV/SVPV.py -example`
* test the gui is working  
`python ./SVPV/SVPV.py -example -gui`
* all done!  

**Note:** SAMtools and BCFtools must be executable by typing 'samtools' and 'bcftools' into the terminal.

###Usage
Running in GUI mode allows users to select and view individual structural variant calls on some subset of the supplied samples. Running in batch mode (i.e. not GUI mode) will generates plots for each call with the suplied set of samples, matching the supplied filter arguments.

**example:**  
```
python SVPV.py -gui -o ./example/output/ -vcf delly:./example/delly.vcf -alt_vcf cnvnator:./example/cnvnator.vcf -manifest ./example/example.manifest -ref_gene ./example/hg38.refgene.partial.txt -ref_vcf ./example/1000G.vcf
```

|Run args:            | Description                                                               | Notes    |
|---------------------|---------------------------------------------------------------------------|----------|
|-vcf<sup>1</sup>     | Primary structural variant prediction VCF/BCF file                        | required |
|-o                   | Output directory                                                          | required |
|-aln                 | Comma separated list of alignment files (indexed BAM/CRAM)                | required <sup>2</sup>
|-samples             | Comma separated list of samples to view, names must be the same as in VCF | required <sup>2</sup>
|-gui                 | run in gui mode                                                           | optional |
|-no_display          | don't attempt to display pdf files in GUI mode                            | optional |
|-alt_vcf<sup>1</sup> | Alternate structural variant prediction vcf/bcf file, <br> called on the same set of samples as primary   | optional
|-ref_vcf<sup>1</sup> | Reference structural variant vcf/bcf file for annotation                  | optional |
|-ref_gene            | Refseq genes regene table file for annotation<sup>3</sup>                 | optional |
|-manifest            | Whitespace delimited file, first column sample names, <br> second column alignment file path. Overrides '-samples' and '-aln' if also given. | optional
|-separate_plots      | Plot each sample separately                                               | optional |

<sup>1</sup>vcfs may be specified by a file (e.g. '-vcf /path/to/file.vcf') or by a name and a file (e.g. '-vcf delly:/path/to/file'). If not specified names will be 'primary', 'alternate' and 'reference' by default.

<sup>2</sup>'-samples' and '-aln' not required if '-manifest' is supplied.

<sup>3</sup>Availble for a variety of reference genomes at [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables)


|Filter args: | Description                                     | Example                          |
------------- |-------------------------------------------------|----------------------------------|
| -max_len    | maximum length of structural variants (bp)      |                                  |
| -min_len    | minimum length of structural variants (bp)      |                                  |
| -af         | Allele frequency threshold                      | -af <0.1                         |
| -gts        | Specify genotypes of given samples              | sample1:0/1,1/1;sample2:1/1      |
| -chrom      | Restrict to comma separated list of chromosomes |                                  |
| -svtype     | Restrict to given SV type (DEL/DUP/CNV/INV)     |                                  |
| -rgi        | Restrict to SVs that intersect refGenes, <br>'-ref_gene' must be supplied          |
| -exonic     | Restrict to SVs that intersect exons of refGenes, <br>'-ref_gene' must be supplied |



|Plot args: | Default | Description                                        |
|-----------|---------|----------------------------------------------------|
|-d         | 1       | force sequencing depth plot on or off              |
|-or        | 1       | force orphaned reads plot on or off                |
|-v         | 1       | force inverted pairs plot on or off                |
|-ss        | 1       | force same strand pairs plot on or off             |
|-hc        | 1       | force hardclipped reads plot on or off             |
|-se        | 0       | force SAM 'secondary alignment' plot on or off     |
|-su        | 0       | force SAM 'supplementary alignment' plot on or off |
|-i         | 1       | force inferred insert size plot on or off          |
|-r         | 1       | force refgenes plot on or off                      |
|-af        | 1       | force allele frequency plot on or off              | 
|-l         | 1       | force plot legend on or off                        |



**another example:**  
```
python SVPV.py -vcf caller1_svs.vcf -samples sample1,sample2,sample3 -aln alignment1.bam,alignment2.bam,alignment3.bam -o /out/directory/ -alt_vcf caller2_svs.vcf -ref_vcf 1000_genomes_svs.vcf -ref_gene hg38.refgene.txt -max_len 100000 -af <0.25 -gts sample1:1/1,0/1;sample3:0/0 -svtype DEL -exonic -ss 0 -se 1
```
