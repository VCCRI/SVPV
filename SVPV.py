# # -*- coding: utf-8 -*-
# """
# author: Jacob Munro, Victor Chang Cardiac Research Institute
# """

from __future__ import print_function
import sys
import os
import re
import config
from svpv.vcf import VCFManager, BCFtools
from svpv.sam import SAMtools
from svpv.refgene import RefgeneManager
from svpv.plot import Plot


def main(argv=sys.argv):
    print("\nStructural Variant Prediction Viewer\n")
    try:
        assert (sys.version_info[0] == 2 and sys.version_info[1] >= 7) or (sys.version_info[0] > 2)
    except AssertionError:
        print("Error: python 2.7+ required. Please update your python installation.")
        exit(1)

    BCFtools.check_installation()
    SAMtools.check_installation()
    if '-example' in argv:
        example(argv)
    else:
        par = Params(argv)
        if not par.run.all:
            par.run.vcf.remove_absent_svs(par.run.samples)

        if par.run.gui:
            import svpv.gui as GUI
            par.filter.gene_list_intersection = False
            GUI.main(par)
        else:
            svs = par.run.vcf.filter_svs(par.filter, as_list=True)
            for sv in svs:
                plot = Plot(sv, par.run.samples, par)
                plot.plot_figure(group=par.plot.grouping)

usage = 'Usage example:\n' \
        'SVPV.py -vcf input_svs.vcf -samples sample1,sample2 -aln alignment1.bam,alignment2.sam -o /out/directory/ \n\n' \
        'Run args:\n' \
        '[required]\n' \
        '-vcf\t\tPrimary structural variant prediction vcf/bcf file.\n' \
        '-aln\t\tComma separated list of alignment files (BAM/CRAM).\n' \
        '-samples\tComma separated list of samples to view,' \
            '\n\t\tNames must be the same as in vcf.\n' \
        '-o\t\tOutput directory.\n' \
        '[optional]\n' \
        '-gui\t\tRun in gui mode.\n' \
        '-no_display\tdon\'t attempt to display pdf files in gui mode \n' \
        '-ref_vcf\tReference structural variant vcf/bcf file for annotation.\n' \
        '-ref_gene\tRefseq genes file for annotation.\n' \
        '-manifest\tWhitespace delimited file with first column sample names and' \
            '\n\t\tsecond column alignment files.' \
            '\n\t\toverrides \'-samples\' and \'-aln\' if also given.\n' \
        'Filter args:\n' \
        '-max_len\tmaximum length of structural variants (bp).\n' \
        '-min_len\tminimum length of structural variants (bp).\n' \
        '-af\t\tAllele frequency threshold, \n' \
            '\t\teg \'-af <0.1\' for SV with frequency less than 0.1.\n' \
        '-gts\t\tSpecify genotypes of given samples.' \
            '\n\t\teg sample1:0/1,1/1;sample3:1/1\n' \
        '-chrom\t\tRestrict to comma separated list of chromosomes.\n' \
        '-svtype\t\tRestrict to given SV type (DEL/DUP/CNV/INV).\n' \
        '-rgi\t\tRestrict to SVs that intersect refGenes, \'-ref_gene\' must be supplied.\n' \
        '-exonic\t\tRestrict to SVs that intersect exons of refGenes,' \
        '\n\t\t\'-ref_gene\' must be supplied.\n' \
        'Plot args:\n' \
        '-d\t0/[1]\tforce sequencing depth plot on or off.\n' \
        '-or\t0/[1]\tforce orphaned reads plot on or off.\n' \
        '-v\t0/[1]\tforce inverted pairs plot on or off.\n' \
        '-ss\t0/[1]\tforce same strand pairs plot on or off.\n' \
        '-cl\t0/[1]\tforce clipped reads plot on or off.\n' \
        '-se\t[0]/1\tforce SAM \'secondary alignment\' plot on or off.\n' \
        '-su\t[0]/1\tforce SAM \'supplementary alignment\' plot on or off.\n' \
        '-i\t0/[1]\tforce inferred insert size plot on or off.\n' \
        '-r\t0/[1]\tforce refgenes plot on or off.\n' \
        '-af\t0/[1]\tforce allele frequency plot on or off.\n' \
        '-l\t0/[1]\tforce plot legend on or off.\n' \
        '\n\nNote: additional parameters can be adjusted by editing \'config.py\'.\n'


def check_file_exists(path):
    if not os.path.isfile(path):
        print("Error: file does not exist!\n'%s'\n" % path)
        exit(1)


class Params:
    def __init__(self, args):
        self.run = RunParams()
        self.filter = FilterParams(self)
        self.plot = PlotParams()

        for i, a in enumerate(args):
            if a[0] == '-':
                # set run parameters
                if a in RunParams.valid:
                    if a == '-vcf':
                        self.run.set_vcfs(args[i + 1])
                    elif a == '-aln':
                        self.run.bams = args[i + 1].split(',')
                    elif a == '-all':
                        self.run.all = True
                    elif a == '-samples':
                        self.run.samples = args[i + 1].split(',')
                    elif a == '-manifest':
                        if self.run.samples or self.run.bams:
                            print("samples and alignments provided as command line arguments overriden by manifest file\n")
                        self.run.read_samples_file(args[i + 1])
                    elif a == '-o':
                        self.run.out_dir = args[i + 1]
                    elif a == '-gui':
                        self.run.gui = True
                    elif a == '-no_display':
                        self.run.display = False
                    elif a == '-ref_gene':
                        self.run.ref_genes = RefgeneManager(args[i + 1])
                        self.filter.ref_genes = self.run.ref_genes
                    elif a == '-ref_vcf':
                        if ':' in args[i + 1]:
                            check_file_exists(args[i + 1].split(':')[1])
                            self.run.ref_vcf = VCFManager(args[i + 1].split(':')[1], name=args[i + 1].split(':')[0], db_mode=True)
                        else:
                            check_file_exists(args[i + 1])
                            self.run.ref_vcf = VCFManager(args[i + 1], name='reference', db_mode=True)

                # set filter parameters
                elif a in FilterParams.valid:
                    # set max and min length filters
                    if a == '-max_len':
                        try:
                            self.filter.max_len = int(args[i + 1])
                        except ValueError:
                            print("invalid max length:" + args[i + 1])
                            exit(1)
                    elif a == '-min_len':
                        try:
                            self.filter.min_len = int(args[i + 1])
                        except ValueError:
                            print("invalid min length:" + args[i + 1])
                            exit(1)
                    # set allele frequency for filtering
                    elif a == '-af':
                        if '>' in args[i + 1]:
                            self.filter.AF_thresh_is_LT = False
                        try:
                            self.filter.AF_thresh = float(re.sub('[<>]', '', args[i + 1]))
                        except ValueError:
                            print(usage)
                            print("invalid allele frequency threshold: -af " + args[i+1])
                            exit(1)
                    # switch for refgene intersections
                    elif a == '-svtype':
                        if args[i+1].upper() in ('DEL', 'DUP', 'INV', 'CNV'):
                            self.filter.svtype = args[i+1].upper()
                        else:
                            print('invalid svtype %s' % args[i+1])
                            exit(1)
                    elif a == '-rgi':
                        self.filter.RG_intersection = True
                    elif a == '-exonic':
                        self.filter.exonic = True
                    # list of genes reported SVs must intersect
                    elif a == '-gene_list':
                        # read in newline/whitespace delimited list of genes
                        self.filter.gene_list = []
                        for line in open(args[i + 1]):
                            for word in line.split():
                                self.filter.gene_list.append(word.strip().upper())
                        self.filter.gene_list_intersection = True
                    elif a == '-gts':
                        # specify genotypes of given samples in form: sample1:0/1,1/1;sample3:1/1
                        self.filter.sampleGTs = {}
                        for sample in args[i + 1].split(';'):
                            filter.sampleGTs[sample.split(':')[0]] = sample.split(':')[1].split(',')
                    elif a == '-chrom':
                        self.filter.chrom = args[i + 1]

                elif a in PlotParams.valid:
                    if a == '-d':
                        if (args[i + 1]) == '0':
                            self.plot.depth = False
                        elif (args[i + 1]) == '1':
                            self.plot.depth = True
                    elif a == '-or':
                        if (args[i + 1]) == '0':
                            self.plot.orphaned = False
                        elif (args[i + 1]) == '1':
                            self.plot.orphaned = True
                    elif a == '-v':
                        if (args[i + 1]) == '0':
                            self.plot.inverted = False
                        elif (args[i + 1]) == '1':
                            self.plot.inverted = True
                    elif a == '-ss':
                        if (args[i + 1]) == '0':
                            self.plot.samestrand = False
                        elif (args[i + 1]) == '1':
                            self.plot.samestrand = True
                    elif a == '-se':
                        if (args[i + 1]) == '0':
                            self.plot.secondary = False
                        elif (args[i + 1]) == '1':
                            self.plot.secondary = True
                    elif a == '-su':
                        if (args[i + 1]) == '0':
                            self.plot.supplementary = False
                        if (args[i + 1]) == '1':
                            self.plot.supplementary = True
                    elif a == '-dm':
                        if (args[i + 1]) == '0':
                            self.plot.diff_mol = False
                        if (args[i + 1]) == '1':
                            self.plot.diff_mol = True
                    elif a == '-cl':
                        if (args[i + 1]) == '0':
                            self.plot.clipped = False
                        elif (args[i + 1]) == '1':
                            self.plot.clipped = True
                    elif a == '-i':
                        if (args[i + 1]) == '0':
                            self.plot.ins = False
                        elif (args[i + 1]) == '1':
                            self.plot.ins = True
                    elif a == '-r':
                        if (args[i + 1]) == '0':
                            self.plot.refgene = False
                        elif (args[i + 1]) == '1':
                            self.plot.refgene = True
                    elif a == '-af':
                        if (args[i + 1]) == '0':
                            self.plot.sv_af = False
                        elif (args[i + 1]) == '1':
                            self.plot.sv_af = True
                    elif a == '-l':
                        if (args[i + 1]) == '0':
                            self.plot.legend = False
                        elif (args[i + 1]) == '1':
                            self.plot.legend = True
                    elif a == '-separate_plots':
                        self.plot.grouping = 1
                else:
                    print("unrecognised argument: " + a)
                    exit(1)
        self.run.check()


# class to store run parameters
class RunParams:
    valid = ('-vcf','-aln', '-samples', '-manifest', '-o', '-gui', '-ref_gene', '-ref_vcf', '-no_display')

    def __init__(self):
        # path to vcf
        self.vcf = None
        # set of alternate sv callsets to visualise against
        self.alt_vcfs = []
        # list of bams
        self.bams = []
        # list of samples
        self.samples = []
        # directory to write data to
        self.out_dir = None
        # switch for gui mode
        self.gui = False
        # refgene manager
        self.ref_genes = None
        # vcf for including population frequencies
        self.ref_vcf = None
        self.all = False

        # get configurations
        # include defaults in case they are accidentally deleted
        try:
            self.display = config.display.split()
        except AttributeError:
            self.display = ['display']
        try:
            self.is_len = config.is_len
        except AttributeError:
            self.is_len = 500
        try:
            self.expansion = config.expansion
        except AttributeError:
            self.expansion = 1
        try:
            self.num_bins = config.num_bins
        except AttributeError:
            self.num_bins = 100

    def get_bams(self, samples):
        bams = []
        for s in samples:
            bams.append(self.bams[self.samples.index(s)])
        return bams

    def read_samples_file(self, filepath):
        check_file_exists(filepath)
        for line in open(filepath):
            if len(line.split()) > 2:
                print("Error: %d fields detected in manifest, explected 2.\n" % len(line.split()))
                exit(1)
            elif len(line.split()) < 2:
                continue
            self.samples.append(line.split()[0].strip())
            self.bams.append(line.split()[1].strip())
            check_file_exists(self.bams[-1])

    # set up the input vcfs (comma separated list, names included with colons name:file or file)
    def set_vcfs(self, vcfs_arg):
        for sv_vcf in vcfs_arg.split(','):
            if ':' in sv_vcf:
                check_file_exists(sv_vcf.split(':')[1])
                vcf = VCFManager(sv_vcf.split(':')[1], name=sv_vcf.split(':')[0])
            else:
                check_file_exists(sv_vcf)
                vcf = VCFManager(sv_vcf)
            if self.vcf is None:
                self.vcf = vcf
            else:
                self.alt_vcfs.append(vcf)

    def check(self):
        if not self.vcf:
            print(usage)
            print("Error: please specify a VCF file")
            exit(1)
        if not self.out_dir:
            print(usage)
            print("Error: please specify out directory")
            exit(1)

        if not self.samples:
            print(usage)
            print("Error: please specify samples to visualise")
            exit(1)
        if not self.bams:
            print(usage)
            print("Error: please specify BAM/SAM files")
            exit(1)
        if not len(self.bams) == len(self.samples):
            print(usage)
            print("Error:\nRequires same number of samples and alignments")
            exit(1)
        for b in self.bams:
            check_file_exists(b)
        delete = []
        for i, s in enumerate(self.samples):
            if s not in self.vcf.samples:
                print("Sample ID not found in VCF: %s - removing from list" % s)
                delete.append(i)
        for i in sorted(delete, reverse=True):
            del self.samples[i]
            del self.bams[i]
        #self.plot_par.check()


# class to store parameters for filtering SVs from VCF
class FilterParams:
    valid = ('-max_len', '-min_len', '-af', '-rgi', '-gene_list', '-gts', '-chrom', '-exonic', '-svtype')

    def __init__(self, parent):
        self.parent = parent
        # threshold for for filtering AF
        self.AF_thresh = None
        # Allele Frequency threshold is less than
        self.AF_thresh_is_LT = True
        # specific chromosome/molecule/contig
        self.chrom = None
        # Dict of list of accepted gentypes for each sample for filtering
        # if sample is not in dict it is not filtered
        self.sample_GTs = {}
        # DEL/DUP/CNV/INV
        self.svtype = None
        # path to genes list file
        self.gene_list = None
        # switch for filtering by gene list
        self.gene_list_intersection = False
        # intersection with refgenes
        self.RG_intersection = False
        # filter SVs by length
        self.min_len = None
        self.max_len = None
        # pointer to refgenes
        self.ref_genes = None
        # filter for SVs that intersect exons only
        self.exonic = False


# class to store parameters for what to show in R plots
class PlotParams:
    valid = ('-d', '-or', '-v', '-ss', '-se', '-su', '-cl', '-i', '-r', '-af', '-l', '-gc', '-dm', '-separate_plots')

    def __init__(self):
        self.gc = False
        self.depth = True
        self.orphaned = True
        self.inverted = True
        self.samestrand = True
        self.secondary = False
        self.supplementary = False
        self.clipped = True
        self.ins = True
        self.refgene = True
        self.sv_af = True
        self.legend = True
        self.diff_mol = True
        self.grouping = 8

    # command line arguments for calling Rscipt
    def get_R_args(self):
        args = []
        if self.depth:
            args.append("-d")
        if self.orphaned:
            args.append("-or")
        if self.inverted:
            args.append("-v")
        if self.samestrand:
            args.append("-ss")
        if self.secondary:
            args.append("-se")
        if self.supplementary:
            args.append("-su")
        if self.clipped:
            args.append("-cl")
        if self.ins:
            args.append("-i")
        if self.refgene:
            args.append("-r")
        if self.sv_af:
            args.append("-af")
        if self.legend:
            args.append("-l")
        if self.gc:
            args.append("-gc")
        if self.diff_mol:
            args.append("-dm")
        return args


def example(argv):
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'example')
    if not os.path.exists(path):
        print('Error: example directory not found')
        exit(1)

    samples = {}
    samples['NA12877_S1'] = os.path.join(path, 'NA12877_S1.partial.bam')
    samples['NA12878_S1'] = os.path.join(path, 'NA12878_S1.partial.bam')
    samples['NA12884_S1'] = os.path.join(path, 'NA12884_S1.partial.bam')
    for s in samples:
        if not os.path.isfile(samples[s]):
            print('Error: example file not found.\n%s' % samples[s])
            exit(1)

    run_argv = []
    for a in argv:
        if not a == '-example':
            run_argv.append(a)
    run_argv.append('-samples')
    run_argv.append(','.join(samples.keys()))
    run_argv.append('-aln')
    run_argv.append(','.join(samples.values()))
    run_argv.append('-vcf')
    run_argv.append('delly:' + os.path.join(path, 'delly.vcf') + ',cnvnator:' + os.path.join(path, 'cnvnator.vcf'))
    run_argv.append('-ref_vcf')
    run_argv.append('1000G:' + os.path.join(path, '1000G.vcf'))
    run_argv.append('-ref_gene')
    run_argv.append(os.path.join(path, 'hg38.refgene.partial.txt'))
    run_argv.append('-o')
    run_argv.append(os.path.join(path, 'output'))
    main(argv=run_argv)
    if '-gui' not in argv:
        print('\nSuccess!\n')

if __name__ == "__main__":
    main()
