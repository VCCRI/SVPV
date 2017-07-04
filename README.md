Structural Variant Prediction Viewer  
------------------------------------
SVPV enables visualisation of predicted structural variant regions in paired-end whole genome sequencing alignments, and
allows comparison of calls from differenct structural variant prediction algorithms. Statistics related to structural
variants are presented in a form that allows users to visually identify false postive calls. Input is a set alignment
files (SAM/BAM/CRAM format) along with a set of structural variant predictions on these alignments (VCF files). Output
is a set of pdf files of structural variant regions. Please see the [wiki](https://github.com/VCCRI/SVPV/wiki/Home) for
examples of SVPV plots of different structural variant categories.

SVPV supports VCF structural variant types deletion (DEL), duplication (DUP), copy number variation (CNV), inversion (INV),
insertion (INS) and breakend ('BND'). Delly2-style translocations (TRA) are also supported.

### Requirements
**Command Line Mode**
* Python 2.7.+ and NumPy
* R v3.+
* [SAMtools and BCFtools](https://github.com/samtools) (version 1.3)
* Linux environment, or access to linux via ssh

**Note:** SAMtools and BCFtools must be executable by typing 'samtools' and 'bcftools' into the terminal.
  
**GUI Mode**
* All command line mode requirements, and:
* X11 if running over ssh
* python 2.7 tkinter
* Recommended: [GraphicsMagick](http://www.graphicsmagick.org/README.html) or ImageMagick ('display')
  * GraphicsMagick build w/ --enable-magick-compat
  * Note: any other X11 capable pdf viewer specified by '-disp' will work
  * alternatively users can navigate to plot directory and open the file with a local pdf viewer

### Installation
* Navigate to desired install directory and clone this repository.

  `git clone https://github.com/VCCRI/SVPV.git SVPV`
* Ensure that requirements are met. For convenience shell scripts are provided for Ubuntu and CentOS.

  `sudo sh ./SVPV/set_up/Ubuntu_set_up.sh`
* Test that SVPV is working. If you get any error messages at this point it is likely that some requirements aren't met.

  `python ./SVPV/SVPV -example`
* Test the gui is working

  `python ./SVPV/SVPV -example -gui`
* all done!

### Non-linux users
* The easiest way to get SVPV running on your Windows or Mac is to run a virtual machine in software such as
 [Oracle VM Virtual Box](https://www.virtualbox.org/).
* You can download an Ubuntu 16.04 image at [osboxes.org](http://www.osboxes.org/ubuntu/)
* After your Ubuntu image is running follow the installation instructions above


### Usage
Running in GUI mode allows users to select and view individual structural variant calls on some subset of the supplied
samples. Running in batch mode (i.e. not GUI mode) will generates plots for each call with the suplied set of samples,
matching the supplied filter arguments.

|Run args:            | Description                                                                | Notes    |
|---------------------|----------------------------------------------------------------------------|----------|
|-vcf<sup>1</sup>     | Comma separated list of structural variant prediction VCF/BCF files        | required |
|-o                   | Output directory                                                           | required |
|-aln                 | Comma separated list of alignment files (indexed BAM/CRAM)                 | required <sup>2</sup>
|-samples             | Comma separated list of samples to view, names must be the same as in VCF  | required <sup>2</sup>
|-gui                 | run in gui mode                                                            | optional |
|-ref_vcf<sup>1</sup> | Reference structural variant vcf/bcf file for annotation                   | optional |
|-ref_gene            | Refseq genes regene table file for annotation<sup>3</sup>                  | optional |
|-manifest            | Whitespace delimited file, first column sample names, <br> second column alignment file path. Overrides '-samples' and '-aln' if also given.                                                     | optional |
|-ped                 | Tab separated pedigree file (GUI only)                                     | optional |
|-separate_plots      | Plot each sample separately                                                | optional |
|-l_svs               | show SVs extending beyond the current plot area.                           | optional |
|-disp                | PDF viewer command. GUI mode only. Default: "display"                      | optional |
|-rd_len              | sequencing read length, optimises window size. Default: 100                | optional |
|-exp                 | window expansion, proportion of SV len added to each side. Default: 1      | optional |
|-bkpt_win            | breakpoint window, number of read lengths to set windows around breakpoints <br> Default:5                                                                                     | optional |
|-n_bins              | target number of bins for plot window. Default: 100                        | optional |



<sup>1</sup>vcfs may be specified by a file (e.g. '-vcf /path/to/file.vcf') or by a name and a file
(e.g. '-vcf delly:/path/to/file'). If not specified names will be 'vcf 1', 'vcf 2', etc and 'reference' by default.

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



|Plot args: | Default | Description                                             |
|-----------|---------|---------------------------------------------------------|
|-d         | 1       | force sequencing depth plot on or off                   |
|-or        | 1       | force orphaned reads plot on or off                     |
|-v         | 1       | force inverted pairs plot on or off                     |
|-ss        | 1       | force same strand pairs plot on or off                  |
|-cl        | 1       | force clipped reads plot on or off                      |
|-se        | 0       | force SAM 'secondary alignment' plot on or off          |
|-su        | 0       | force SAM 'supplementary alignment' plot on or off      |
|-dm        | 0       | force mate different molecule alignment plot on or off  |
|-i         | 1       | force inferred insert size plot on or off               |
|-r         | 1       | force refgenes plot on or off                           |
|-af        | 1       | force allele frequency plot on or off                   |
|-l         | 1       | force plot legend on or off                             |

**Usage example:**
```
python SVPV -gui -o ./example/output/ -vcf delly:./example/delly.vcf,cnvnator:./example/cnvnator.vcf -manifest
/example/example.manifest -ref_gene ./example/hg38.refgene.partial.txt -ref_vcf 1000G:./example/1000G.vcf
```
**Advance usage example:**
```
python SVPV -vcf caller1_svs.vcf,caller2_svs.vcf -samples sample1,sample2,sample3 -aln s1.bam,s2.bam,s3.bam
-o /out/directory/ -ref_vcf 1000_genomes_svs.vcf -ref_gene hg38.refgene.txt -max_len 100000 -af <0.25 -gts
sample1:1/1,0/1;sample3:0/0 -svtype DEL -exonic -ss 0 -se 1
```

###  VCF Field Requirements:

SV Type         | Required VCF Fields
----------------|------------------------------------------------|
All Types       | CHROM, POS, SVTYPE<sup>1</sup>, GT<sup>2</sup> |
DEL/DUP/CNV/INV | SVLEN or END                                   |
INS             | SVLEN or INSLEN\*                              |
BND             | ALT, ID, MATEID/PAIRID/EVENTID                 |
TRA\*           | CHR2, END                                      |


<sup>1</sup>If SVTYPE is not found then ALT is parsed for symbolic alternate alleles.
These should match one of DEL, DUP, CNV, INS, INV, BND or TRA or the call will be ignored.

<sup>2</sup>For reference VCF only, GT is not required if AF is present.

*Included for compatibility with Delly2

Please see the [VCF specifications](http://samtools.github.io/hts-specs/VCFv4.3.pdf) for further details.

### SVPV Plot Window Sizing and Types
* Please see the [wiki](https://github.com/VCCRI/SVPV/wiki/SVPV-Plot-Windows)
