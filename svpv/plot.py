# # -*- coding: utf-8 -*-
# """
# author: Jacob Munro, Victor Chang Cardiac Research Institute
# """
from __future__ import print_function
from __future__ import division
import os
import subprocess
from hashlib import sha1
import copy
from .sam import SamStats, SAMtools
from .vcf import SV
from .refgene import RefGeneEntry


class Plot:
    svpv_r = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'svpv.r')

    def __init__(self, sv, samples, par):
        self.par = par
        self.samples = samples
        self.sv = sv
        self.region_bins = None
        self.bkpt_bins = None

        # half breakpoint window
        h_bkpt_wind = (par.run.bkpt_win * par.run.rd_len) // 2

        # show depth over whole region but zoom in on breakpoints if necessary
        if sv.svtype in ('DEL', 'DUP', 'CNV', 'INV', 'CUSTOM'):
            start = sv.pos - par.run.expansion * (sv.end - sv.pos + 1)
            end = sv.end + par.run.expansion * (sv.end - sv.pos + 1)
            self.region_bins = Bins(sv.chrom, start, end, ideal_num_bins=par.run.num_bins)
            if (end - start) > self.par.run.bkpt_win * self.par.run.rd_len:
                self.bkpt_bins = (Bins(sv.chrom, sv.pos - h_bkpt_wind, sv.pos + h_bkpt_wind,
                                       ideal_num_bins=par.run.num_bins//2),
                                  Bins(sv.chrom, sv.end - h_bkpt_wind, sv.end + h_bkpt_wind,
                                       ideal_num_bins=par.run.num_bins//2))
                self.sam_stats = SamStats.get_sam_stats(par.run.get_bams(samples), self.bkpt_bins,
                                                   depth_bins=self.region_bins)
            elif self.region_bins.length() <  self.par.run.bkpt_win * self.par.run.rd_len:
                mid = (sv.pos + sv.end) // 2
                self.region_bins = Bins(sv.chrom, mid - h_bkpt_wind, mid + h_bkpt_wind, ideal_num_bins=par.run.num_bins)
                self.sam_stats = SamStats.get_sam_stats(par.run.get_bams(samples), [self.region_bins])
            else:
                self.sam_stats = SamStats.get_sam_stats(par.run.get_bams(samples), [self.region_bins])

        # single breakpoint
        elif sv.svtype == 'INS':
            self.region_bins = Bins(sv.chrom, sv.pos - h_bkpt_wind, sv.pos + h_bkpt_wind,
                                    ideal_num_bins=par.run.num_bins)
            self.sam_stats = SamStats.get_sam_stats(par.run.get_bams(samples), [self.region_bins])

        # do not show depth region, just stats at pair of breakpoints
        elif sv.svtype in ('BND', 'TRA'):
            if sv.svtype == 'BND':
                chr1, pos1 = sv.BND_Event.loci[0]
                if len(sv.BND_Event.loci) > 1:
                    chr2, pos2 = sv.BND_Event.loci[1]
                else:
                    chr2, pos2 = None, None
            # delly TRA exception
            elif sv.svtype == 'TRA':
                chr1, pos1 = sv.chrom, sv.pos
                chr2, pos2 = sv.chr2, sv.chr2_pos
            if chr1 == chr2 and abs(pos2-pos1) < 2*par.run.rd_len:
                self.region_bins = Bins(chr1, pos1 - h_bkpt_wind, pos2 + h_bkpt_wind, ideal_num_bins=par.run.num_bins)
                self.sam_stats = SamStats.get_sam_stats(par.run.get_bams(samples), [self.region_bins])
            elif chr2 is not None:
                self.bkpt_bins = (Bins(chr1, pos1 - h_bkpt_wind, pos1 + h_bkpt_wind, ideal_num_bins=par.run.num_bins//2),
                                  Bins(chr2, pos2 - h_bkpt_wind, pos2 + h_bkpt_wind, ideal_num_bins=par.run.num_bins//2))
                self.sam_stats = SamStats.get_sam_stats(par.run.get_bams(samples), self.bkpt_bins)
            else:
                self.region_bins = Bins(chr1, pos1 - h_bkpt_wind, pos1 + h_bkpt_wind, ideal_num_bins=par.run.num_bins)
                self.sam_stats = SamStats.get_sam_stats(par.run.get_bams(samples), [self.region_bins])

        else:
            raise ValueError('unsupported svtype: {}'.format(self.sv.svtype))
        self.print_data()

    def print_data(self):
        # create directories
        self.dirs = self.create_dirs(self.par.run.out_dir)

        # print sample data to file
        for i, s in enumerate(self.samples):
            self.sam_stats[i].print_stats(self.dirs[s])

        if self.par.run.fa:
            if self.region_bins:
                self.region_bins.print_gc(self.par.run.fa, open(os.path.join(self.dirs['pos'], 'region_gc.tsv'), 'wt'))
            else:
                for bin in self.bkpt_bins:
                    bin.print_gc(self.par.run.fa, open(os.path.join(self.dirs['pos'],
                                                                    '{}.{}.gc.tsv'.format(bin.chrom, bin.start)), 'wt'))

        # plot attributes for use in R
        plot_attr = open(os.path.join(self.dirs['pos'], 'plot_attr.tsv'), 'wt')
        plot_attr.write('\t'.join(('ver', 'region', 'r_bin_size', 'r_bin_num', 'loci', 'l_bin_size', 'l_bin_num')) + '\n')
        plot_attr.write('SVPV v{}\t'.format(self.par.ver))
        if self.region_bins:
            plot_attr.write('{}\t{}\t{}\t'.format(self.region_bins.region, self.region_bins.size, self.region_bins.num))
        else:
            plot_attr.write('NA\tNA\tNA\t')
        if self.bkpt_bins:
            for bin in self.bkpt_bins:
                plot_attr.write('{},'.format(bin.region))
            plot_attr.write('\t{}\t{}\n'.format(self.bkpt_bins[0].size, self.bkpt_bins[0].num))
        else:
            plot_attr.write('NA\tNA\tNA\n')
        plot_attr.close()

        # extract query regions
        if self.region_bins:
            queries = [self.region_bins.get_region_tuple()]
        else:
            queries = []
            for bin in self.bkpt_bins:
                queries.append(bin.get_region_tuple())

        # gene annotation
        if self.par.run.ref_genes:
            genes = []
            for region in queries:
                gs =  self.par.run.ref_genes.get_entries_in_range(*region)
                for g in gs:
                    if g not in genes:
                        genes.append(g)
            if genes:
                RefGeneEntry.print_entries(genes, open(os.path.join(self.dirs['pos'], 'refgene.tsv'), 'w'))

        # sample-wise SV annotation
        for i, s in enumerate(self.samples):
            sv_file = open(os.path.join(self.dirs[s], 'svs.tsv'), 'w')
            svs = []
            for region in queries:
                _svs_ = self.par.run.vcf.get_svs_in_range(*region, sample=s, lrg_svs=self.par.plot.l_svs)
                for sv in _svs_:
                    if sv not in svs:
                        svs.append(sv)
            SV.print_SVs_header(sv_file, sample_index=self.par.run.vcf.get_sample_index(s))
            SV.print_SVs(svs, sv_file, self.par.run.vcf.name, sample_index=self.par.run.vcf.get_sample_index(s))
            for vcf in self.par.run.alt_vcfs:
                svs = []
                for region in queries:
                    _svs_ = vcf.get_svs_in_range(*region, sample=s, lrg_svs=self.par.plot.l_svs)
                    for sv in _svs_:
                        if sv not in svs:
                            svs.append(sv)
                if svs:
                    SV.print_SVs(svs, sv_file, vcf.name, sample_index=vcf.get_sample_index(s))
            sv_file.close()

        # batch-wise SV annotation
        svs_file = open(os.path.join(self.dirs['pos'], 'SV_AF.tsv'), 'w')
        SV.print_SVs_header(svs_file)
        # primary vcf
        svs = []
        for region in queries:
            _svs_ = self.par.run.vcf.get_svs_in_range(*region, lrg_svs=self.par.plot.l_svs)
            for sv in _svs_:
                if sv not in svs:
                    svs.append(sv)
        if svs:
            SV.print_SVs(svs, svs_file, self.par.run.vcf.name)
        if self.sv.svtype == 'CUSTOM':
            SV.print_SVs([self.sv], svs_file, 'CUSTOM')
        # ref vcf
        if self.par.run.ref_vcf:
            svs = []
            for region in queries:
                _svs_ = self.par.run.ref_vcf.get_svs_in_range(*region, lrg_svs=self.par.plot.l_svs)
                for sv in _svs_:
                    if sv not in svs:
                        svs.append(sv)
            if svs:
                SV.print_SVs(svs, svs_file, self.par.run.ref_vcf.name)
        # alt vcfs
        for vcf in self.par.run.alt_vcfs:
            svs = []
            for region in queries:
                _svs_ = vcf.get_svs_in_range(*region, lrg_svs=self.par.plot.l_svs)
                for sv in _svs_:
                    if sv not in svs:
                        svs.append(sv)
            if svs:
                SV.print_SVs(svs, svs_file, vcf.name)
        svs_file.close()

    def plot_figure(self, group=8, display=False):
        # split into groups of 8 or less so don't go over R layout limit
        out = ''
        current_samples = self.samples[0:group]
        next_samples = self.samples[group:]
        while True:
            if group == 1:
                id = current_samples[0]
            else:
                id = sha1(''.join(current_samples).encode('utf-8')).hexdigest()[0:10]
            out = os.path.join(self.dirs['pos'], '{}.{}.{}.{}.pdf'.format(self.sv.chrom, self.sv.pos, self.sv.svtype,id))
            cmd = ['Rscript', Plot.svpv_r, ','.join(current_samples), os.path.join(self.dirs['pos'], ''), out]
            cmd.append('"{} at {}:{}"'.format(self.sv.svtype, self.sv.chrom, self.sv.pos))
            cmd.extend(self.par.plot.get_R_args())
            print(' '.join(cmd) + '\n')
            try:
                subprocess.check_call(cmd)
            except OSError:
                print('Rscript failed. Are you sure it is installed?')
                exit(1)

            if display:
                cmd = [display]
                cmd.append(out)
                print(' '.join(cmd) + '\n')
                try:
                    subprocess.check_call(cmd)
                except OSError:
                    print('Error: could not run %s. Are you sure it is installed?' % ' '.join(display))
                    exit(1)
            else:
                print("created %s\n" % out)
            current_samples = next_samples[0:group]
            next_samples = next_samples[group:]
            if not current_samples:
                break
        return out

    def create_dirs(self, outdir):
        dirs = {}
        dirs['root'] = outdir
        if not os.path.exists(dirs['root']):
            os.makedirs(dirs['root'])
        dirs['svtype'] = os.path.join(dirs['root'], self.sv.svtype)
        if not os.path.exists(dirs['svtype']):
            os.mkdir(dirs['svtype'])
        dirs['pos'] = os.path.join(dirs['svtype'], '%s_%d' % (self.sv.chrom, self.sv.pos))
        if not os.path.exists(dirs['pos']):
            os.mkdir(dirs['pos'])
        for s in self.samples:
            dirs[s] = os.path.join(dirs['pos'], s)
            if not os.path.exists(dirs[s]):
                os.mkdir(dirs[s])
        return dirs

class Bins:
    def __init__(self, chrom, start, end, ideal_num_bins=100):

        # aim for ideal_num_bins, but bins need to be uniformally distributed and of equal size
        # smallest bins size is 1bp, so for regions < num_bins bp there will be less than num_bins bins
        self.chrom = chrom
        self.start = max(0, start)
        self.size = (end - self.start + 1) // ideal_num_bins
        self.size += not (self.size) * 1
        self.num = (end - self.start + 1) // self.size
        self.end = self.start + self.num * self.size - 1
        self.region = chrom + ':' + str(self.start) + '-' + str(self.end)

    def print_gc(self, fasta, file):
        vals =[]
        for i in range(self.num):
            vals.append(str(SAMtools.get_GC(fasta, '{}:{}-{}'.format(self.chrom, self.start + i*self.size,
                                                                                self.start + (i+1)*self.size))))
        file.write('\t'.join(vals) + '\n')
        file.close()

    def length(self):
        return self.end - self.start + 1

    def get_region_tuple(self):
        return self.chrom, self.start, self.end

    def get_bin_coverage(self, start, end):
        if start > self.end or end < self.start:
            return None

        if start <= self.start:
            first = 0
            if end >= self.start + self.size:
               first_bp = self.size
            else:
               first_bp = (end - self.start + 1)
        else:
            first = (start- self.start) // self.size
            first_bp = self.size - ((start - self.start + 1) % self.size)

        if end >= self.end:
            last = self.num - 1
            if start <= self.end - self.size:
                last_bp = self.size
            else:
                last_bp = start- (self.end - self.size)
        else:
            last = (end - self.start) // self.size
            last_bp = (end - self.start + 1) % self.size

        if first < 0 or first >= self.num or last < 0 or last >= self.num:
            return None

        if first == last:
            last_bp = 0
        return (first, first_bp), (last, last_bp)
