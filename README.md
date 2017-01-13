Structural Variant Prediction Viewer  
------------------------------------
View predicted structural variant regions from paired-end whole genome sequencing alignments and compare calls from
differenct structural variant prediction algorithms. Statistics related to structural variants are presented in a form
that allows users to visually identify false postive calls. Input is a set alignment files (SAM/BAM/CRAM format), and a
set of structural variant predictions on these alignments (VCF files). Output is a set of pdf files of structural
variant regions.

VCF structural variant types deletion (DEL), duplication (DUP), copy number variation (CNV), inversion (INV),
insertion (INS) and breakend ('BND') are supported. Delly-style translocations (TRA) are also supported.

###Requirements
**Basic**  
* Python 2.7.+ and NumPy
* R v3.+
* [SAMtools and BCFtools](https://github.com/samtools) (version 1.3)
* Linux environment, or access to linux via ssh
**Note:** SAMtools and BCFtools must be executable by typing 'samtools' and 'bcftools' into the terminal.
  
**GUI**  
* X11 if running over ssh
* python 2.7 tkinter
* [GraphicsMagick](http://www.graphicsmagick.org/) or ImageMagick ('display')
  * Recommended: GraphicsMagick build w/ --enable-magick-compat
  * or other X11 pdf viewer specified in 'config.py'

###Installation
* Navigate to desired install directory and clone this repository.

  `git clone https://github.com/VCCRI/SVPV.git SVPV`
* Ensure that requirements are met. For convenience shell scripts are provided for Ubuntu and CentOS.

  `sudo sh ./SVPV/set_up/Ubuntu_set_up.sh`
* Test that SVPV is working. If you get any error messages at this point it is likely that some requirements aren't met.

  `python ./SVPV/SVPV -example`
* Test the gui is working

  `python ./SVPV/SVPV -example -gui`
* all done!

###Non-linux users
* The easiest way to get SVPV running on your PC or Mac is to run a virtual machine in software such as
 [Oracle VM Virtual Box](https://www.virtualbox.org/).
* You can download an Ubuntu 16.04 image at [osboxes.org](http://www.osboxes.org/ubuntu/)
* After your Ubuntu image is running follow the installation instructions above


###Usage
Running in GUI mode allows users to select and view individual structural variant calls on some subset of the supplied samples. Running in batch mode (i.e. not GUI mode) will generates plots for each call with the suplied set of samples, matching the supplied filter arguments.

**example:**  
```
python SVPV -gui -o ./example/output/ -vcf delly:./example/delly.vcf -alt_vcf cnvnator:./example/cnvnator.vcf -manifest ./example/example.manifest -ref_gene ./example/hg38.refgene.partial.txt -ref_vcf ./example/1000G.vcf
```
**another example:**
```
python SVPV -vcf caller1_svs.vcf -samples sample1,sample2,sample3 -aln alignment1.bam,alignment2.bam,alignment3.bam -o /out/directory/ -alt_vcf caller2_svs.vcf -ref_vcf 1000_genomes_svs.vcf -ref_gene hg38.refgene.txt -max_len 100000 -af <0.25 -gts sample1:1/1,0/1;sample3:0/0 -svtype DEL -exonic -ss 0 -se 1
```

|Run args:            | Description                                                               | Notes    |
|---------------------|---------------------------------------------------------------------------|----------|
|-vcf<sup>1</sup>     | list of structural variant prediction VCF/BCF files                       | required |
|-o                   | Output directory                                                          | required |
|-aln                 | Comma separated list of alignment files (indexed BAM/CRAM)                | required <sup>2</sup>
|-samples             | Comma separated list of samples to view, names must be the same as in VCF | required <sup>2</sup>
|-gui                 | run in gui mode                                                           | optional |
|-no_display          | don't attempt to display pdf files in GUI mode                            | optional |
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



### Structural Variant VCFs
BCFtools is used to parse vcf/bcf formatted files.
SVPV expects either the 'SVTYPE' info field or symbolic alternative alleles (e.g. '\<DEL\>') to recognise structural variant calls.
VCF entries with neither 'SVTYPE' or symbolic allele on the supported SVtype list (DEL, DUP, CNV, INS\*, INV, BND, TRA\*) will be ignored.
'END' is required for deletion, duplication and CNV type variants and 'ISLEN' for insertions.
SVTYPE='BND': Translocations, inversions and generic breakend types are also supported. These require the either MATEID or EVENTID fields.
Please see the [VCF specifications](http://samtools.github.io/hts-specs/VCFv4.3.pdf) for clarification.

\*For compatipility with Delly, SVTYPE='TRA' is supported, and info fields 'CHR2' and 'INSLEN' are parsed.