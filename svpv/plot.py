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
from .SAM import SamStats
from .VCF import SV
from .refgene import RefGeneEntry


class Plot:
    svpv_r = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'SVPV.r')

    def __init__(self, sv, samples, par, breakpoint_zoom=True):
        self.par = par
        self.sv = sv
        self.start = sv.start - par.run.expansion * (sv.end - sv.start + 1)
        self.end = sv.end + par.run.expansion * (sv.end - sv.start + 1)
        self.samples = samples
        self.dirs = self.create_dirs(par.run.out_dir)

        bin_size = (self.end - self.start) // par.run.num_bins
        if breakpoint_zoom and bin_size // float(par.run.is_len) > 0.25:
            breakpoints = ((sv.start - int(1.5 * par.run.is_len), sv.start + int(1.5 * par.run.is_len)),
                           (sv.end - int(1.5 * par.run.is_len), sv.end + int(1.5 * par.run.is_len)))
            sam_stats = SamStats.get_sam_stats(sv.chrom, self.start, self.end, par.run.get_bams(samples),
                                               par.run.num_bins, breakpoints=breakpoints)
        else:
            sam_stats = SamStats.get_sam_stats(sv.chrom, self.start, self.end, par.run.get_bams(samples),
                                               par.run.num_bins)

        # print sample data to file
        for i, s in enumerate(samples):
            sam_stats[i].print_stats(self.dirs[s])

        if par.run.ref_genes:
            genes = par.run.ref_genes.get_entries_in_range(sv.chrom, self.start, self.end)
            if genes:
                RefGeneEntry.print_entries(genes, open(os.path.join(self.dirs['pos'], 'refgene.tsv'), 'w'))

        self.add_vcf_annotation()

    def add_vcf_annotation(self):
        # sample-wise SV annotation
        for i, s in enumerate(self.samples):
            sv_file = open(os.path.join(self.dirs[s], 'svs.tsv'), 'w')
            svs = self.par.run.vcf.get_svs_in_range(self.sv.chrom, self.start, self.end, sample=s)
            if self.sv not in svs:
                svs.append(self.sv)
            SV.print_SVs_header(sv_file, sample_index=self.par.run.vcf.get_sample_index(s))
            SV.print_SVs(svs, sv_file, self.par.run.vcf.name, sample_index=self.par.run.vcf.get_sample_index(s))
            for vcf in self.par.run.alt_vcfs:
                alt_svs = vcf.get_svs_in_range(self.sv.chrom, self.start, self.end, sample=s)
                if alt_svs:
                    SV.print_SVs(alt_svs, sv_file, vcf.name, sample_index=vcf.get_sample_index(s))
            sv_file.close()

        # batch-wise SV annotation
        svs_file = open(os.path.join(self.dirs['pos'], 'SV_AF.tsv'), 'w')
        SV.print_SVs_header(svs_file)
        batch_svs = self.par.run.vcf.get_svs_in_range(self.sv.chrom, self.start, self.end)
        if batch_svs:
            SV.print_SVs(batch_svs, svs_file, self.par.run.vcf.name)
        if self.par.run.ref_vcf:
            ref_svs = self.par.run.ref_vcf.get_svs_in_range(self.sv.chrom, self.start, self.end)
            if ref_svs:
                SV.print_SVs(ref_svs, svs_file, self.par.run.ref_vcf.name)
        for vcf in self.par.run.alt_vcfs:
            alt_svs = vcf.get_svs_in_range(self.sv.chrom, self.start, self.end)
            if alt_svs:
                SV.print_SVs(alt_svs, svs_file, vcf.name)
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
            out = os.path.join(self.dirs['pos'], '%s.%s.%s.%s.%s.pdf' % (self.sv.chrom, self.sv.start, self.sv.svtype,
                                                                         self.get_length_units(), id))
            cmd = ['Rscript']
            cmd.append(Plot.svpv_r)
            cmd.append(','.join(current_samples))
            cmd.append(os.path.join(self.dirs['pos'], ''))
            cmd.append(out)
            cmd.append('"%s at %s:%d-%d"' % (self.sv.svtype, self.sv.chrom, self.sv.start, self.sv.end))
            cmd.extend(self.par.plot.get_R_args())
            print(' '.join(cmd) + '\n')
            try:
                subprocess.call(cmd)
            except OSError:
                print('Error: failed to run Rscript. Are you sure R is installed?')
                exit(1)

            if display:
                cmd = copy.copy(display)
                cmd.append(out)
                print(' '.join(cmd) + '\n')
                try:
                    subprocess.call(cmd)
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
        dirs['pos'] = os.path.join(dirs['svtype'], '%s_%d-%d' % (self.sv.chrom, self.sv.start, self.sv.end))
        if not os.path.exists(dirs['pos']):
            os.mkdir(dirs['pos'])
        for s in self.samples:
            dirs[s] = os.path.join(dirs['pos'], s)
            if not os.path.exists(dirs[s]):
                os.mkdir(dirs[s])
        return dirs

    def get_length_units(self):
        length = self.sv.end - self.sv.start + 1
        if length < 1e3:
            return '%d_%s' % (length, 'bp')
        elif 1e3 <= length < 1e6:
            return '%d_%s' % (length/1e3, 'kbp')
        elif 1e6 <= length < 1e9:
            return '%d_%s' % (length/1e6, 'Mbp')
        else:
            return '%d_%s' % (length/1e9, 'Gbp')
