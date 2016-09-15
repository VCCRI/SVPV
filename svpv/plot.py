# # -*- coding: utf-8 -*-
# """
# author: Jacob Munro, Victor Chang Cardiac Research Institute
# """

import os
import subprocess
from hashlib import sha1
import copy
from SAM import SamStats
from VCF import SV
from refgene import RefGeneEntry


class Plot:
    svpv_r = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'SVPV.r')

    def __init__(self, sv, samples, par, expansion=1, num_bins=100, is_len=500, breakpoint_zoom=True):
        self.par = par
        self.sv = sv
        self.start = sv.start - expansion * (sv.end - sv.start + 1)
        self.end = sv.end + expansion * (sv.end - sv.start + 1)
        self.samples = samples
        self.dirs = self.create_dirs(par.run.out_dir)

        bin_size = (self.end - self.start)/num_bins
        if breakpoint_zoom and bin_size / float(is_len) > 0.25:
            breakpoints = ((sv.start - int(1.5*is_len), sv.start + int(1.5*is_len)), (sv.end - int(1.5*is_len),
                          sv.end + int(1.5*is_len)))
            sam_stats = SamStats.get_sam_stats(sv.chrom, self.start, self.end, par.run.get_bams(samples), num_bins,
                                               breakpoints=breakpoints)
        else:
            sam_stats = SamStats.get_sam_stats(sv.chrom, self.start, self.end, par.run.get_bams(samples), num_bins)

        # print sample data to file
        for i, s in enumerate(samples):
            sam_stats[i].print_stats(self.dirs[s])

        if par.run.ref_genes:
            genes = par.run.ref_genes.get_entries_in_range(sv.chrom, self.start, self.end)
            if genes:
                RefGeneEntry.print_entries(genes, file(os.path.join(self.dirs['pos'], 'refgene.tsv'), 'w'))

        self.add_vcf_annotation()

    def add_vcf_annotation(self):
        for i, s in enumerate(self.samples):
            sv_file = file(os.path.join(self.dirs[s], 'svs.tsv'), 'w')
            svs = self.par.run.vcf.get_svs_in_range(self.sv.chrom, self.start, self.end, sample=s)
            if self.sv not in svs:
                svs.append(self.sv)
            SV.print_SVs_header(sv_file, sample_index=self.par.run.vcf.get_sample_index(s))
            SV.print_SVs(svs, sv_file, self.par.run.vcf.name, sample_index=self.par.run.vcf.get_sample_index(s))
            if self.par.run.alt_vcf:
                alt_svs = self.par.run.alt_vcf.get_svs_in_range(self.sv.chrom, self.start, self.end, sample=s)
            else:
                alt_svs = None
            if alt_svs:
                SV.print_SVs(alt_svs, sv_file, self.par.run.alt_vcf.name,
                             sample_index=self.par.run.alt_vcf.get_sample_index(s))
            sv_file.close()

        # print annotation data to file
        batch_svs = self.par.run.vcf.get_svs_in_range(self.sv.chrom, self.start, self.end)
        if self.par.run.ref_vcf:
            ref_svs = self.par.run.ref_vcf.get_svs_in_range(self.sv.chrom, self.start, self.end)
        else:
            ref_svs = None
        if self.par.run.alt_vcf:
            alt_svs = self.par.run.alt_vcf.get_svs_in_range(self.sv.chrom, self.start, self.end)
        else:
            alt_svs = None

        svs_file = file(os.path.join(self.dirs['pos'], 'SV_AF.tsv'), 'w')
        SV.print_SVs_header(svs_file)
        if batch_svs:
            SV.print_SVs(batch_svs, svs_file, self.par.run.vcf.name)
        if alt_svs:
            SV.print_SVs(alt_svs, svs_file, self.par.run.alt_vcf.name)
        if ref_svs:
            SV.print_SVs(ref_svs, svs_file, self.par.run.ref_vcf.name)
        svs_file.close()

    def plot_figure(self, display=False):
        # split into groups of 8 or less so don't go over R layout limit
        out = ''
        current_samples = self.samples[0:8]
        next_samples = self.samples[8:]
        while True:
            sha_id = sha1(''.join(current_samples)).hexdigest()[0:10]
            out = os.path.join(self.dirs['pos'], '%s.%s.%s.%s.%s.pdf' % (self.sv.chrom, self.sv.start, self.sv.svtype,
                                                                         self.get_good_length_units(), sha_id))
            cmd = ["Rscript"]
            cmd.append(Plot.svpv_r)
            cmd.append(','.join(current_samples))
            cmd.append(os.path.join(self.dirs['pos'], ''))
            cmd.append(out)
            cmd.append('"%s at %s:%d-%d"' % (self.sv.svtype, self.sv.chrom, self.sv.start, self.sv.end))
            cmd.extend(self.par.plot.get_R_args())
            print ' '.join(cmd) + '\n'

            try:
                subprocess.call(cmd)
            except OSError:
                print 'Error: could not run Rscript. Are you sure it is installed?'
                exit(1)

            if display:
                cmd = copy.copy(display)
                cmd.append(out)
                try:
                    subprocess.call(cmd)
                except OSError:
                    print 'Error: could not run %s. Are you sure it is installed?' % ' '.join(display)
                    exit(1)
                print ' '.join(cmd) + '\n'
            else:
                print "created %s\n" % out
            current_samples = next_samples[0:8]
            next_samples = next_samples[8:]
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

    def get_good_length_units(self):
        length = self.sv.end - self.sv.start + 1
        if length < 1e3:
            return '%d_%s' % (length, 'bp')
        elif 1e3 <= length < 1e6:
            return '%d_%s' % (length/1e3, 'kbp')
        elif 1e6 <= length < 1e9:
            return '%d_%s' % (length/1e6, 'Mbp')
        else:
            return '%d_%s' % (length/1e9, 'Gbp')


