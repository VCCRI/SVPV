# # -*- coding: utf-8 -*-
# """
# author: Jacob Munro, Victor Chang Cardiac Research Institute
# """

import re
import subprocess
import os
from subprocess import PIPE
import numpy as np


class SamEntry():
    # regex for processing cigar string
    cigar_hard_clip = re.compile('^((?P<L>[0-9]+)H)?([0-9]+[^H])+((?P<R>[0-9]+)H)?$')
    cigar_ref_chars = re.compile('[SMDX=N]')
    # sam flags
    paired = np.uint16(1)
    mapped_in_proper_pair = np.uint16(2)
    read_unmapped = np.uint16(4)
    mate_unmapped = np.uint16(8)
    read_reverse = np.uint16(16)
    mate_reverse = np.uint16(32)
    first = np.uint16(64)
    second = np.uint16(128)
    secondary = np.uint16(256)
    fails_QC = np.uint16(512)
    duplicate = np.uint16(1024)
    supplementary = np.uint16(2048)

    def __init__(self, fields):
        self.flag = np.uint16(fields[1])
        self.pos = int(fields[3])
        self.mapQ = int(fields[4])
        self.cigar = fields[5]
        self.tlen = int(fields[8])

    def has_flag(self, flag):
        return self.flag & flag
    def is_rvs(self):
        return (self.flag & SamEntry.read_reverse)

    def has_unmapped_mate(self):
        return self.flag & SamEntry.mate_unmapped

    # not_primary or supplementary alignment
    def is_alt_alignment(self):
        return (self.flag & SamEntry.supplementary) or (self.flag & SamEntry.not_primary)

    def mate_same_strand(self):
        return ((self.flag & SamEntry.read_reverse) and (self.flag & SamEntry.mate_reverse)) or \
               (not(self.flag & SamEntry.read_reverse)) and (not(self.flag & SamEntry.mate_reverse))

    # assume not same strand
    def is_inverted(self):
        if not(self.flag & SamEntry.read_reverse):
            if self.tlen < 0:
                return True
        elif self.tlen > 0:
            return True
        return False

    # return the number of ref bases between the first and last mapped positions
    def get_right_pos(self):
        num = 0
        num_s = ''
        for c in self.cigar:
            if c.isdigit():
                num_s += c
            else:
                if re.match(self.cigar_ref_chars, c):
                    num += int(num_s)
                num_s = ''
        return num + self.pos - 1

    # assume entry and mate are mapped
    def is_discordant(self):
        if (self.flag & SamEntry.read_reverse):
            if(self.flag & SamEntry.mate_reverse):
                '''  <--1 <--2  or  <--2 <--1  '''
                return True
            elif self.tlen > 0:
                '''  <--1 2-->  '''
                return True
            else:
                '''  2--> <--1  '''
                return False
        elif not (self.flag & SamEntry.read_reverse):
            if not (self.flag & SamEntry.mate_reverse):
                '''  1--> 2-->  or  2--> 1-->  '''
                return True
            elif self.tlen < 0:
                '''  <--2 1-->  '''
                return True
            else:
                '''  1--> <--2  '''
                return False
        return False

    def get_hard_clips(self):
        m = re.search(SamEntry.cigar_hard_clip, self.cigar)
        if m.group('R'):
            if m.group('L'):
                return (int(m.group('L')), int(m.group('R')))
            else:
                return (int(0), int(m.group('R')))
        elif m.group('L'):
            return (int(m.group('L')), int(0))
        else:
            return (0,0)


class Bins():
    def __init__(self, start, end, ideal_num_bins=100):
        # aim for ideal_num_bins, but bins need to be uniformally distributed and of equal size
        # smallest bins size is 1bp, so for regions < num_bins bp there will be less than num_bins bins
        self.start = start
        self.size = (end - start + 1) / ideal_num_bins
        self.size += not (self.size) * 1
        self.num = (end - start + 1) / self.size
        self.end = self.start + self.num * self.size - 1

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
            first = (start- self.start) / self.size
            first_bp = self.size - ((start - self.start + 1) % self.size)

        if end >= self.end:
            last = self.num - 1
            if start <= self.end - self.size:
                last_bp = self.size
            else:
                last_bp = start- (self.end - self.size)
        else:
            last = (end - self.start) / self.size
            last_bp = (end - self.start + 1) % self.size

        if first < 0 or first >= self.num or last < 0 or last >= self.num:
            return None

        if first == last:
            last_bp = 0
        return ((first, first_bp),(last,last_bp))


class SamStats():
    depth_cols = ['total', 'mapQltT', 'mapQ0']
    TOTAL = 0
    MAPQLTT = 1
    MAPQ0 = 2
    aln_stats_cols = ['reads', 'orphaned', 'inverted', 'samestrand', 'secondary', 'supplementary', 'hardclipped']
    READS = 0
    ORPHANED = 1
    INVERTED = 2
    SAMESTRAND = 3
    SECONDARY = 4
    SUPPLEMENTARY = 5
    HARDCLIPPED = 6

    def __init__(self, depth_bins, bkpt_bins=None, mapq_t=30, hardclip_t=1):
        # set parameters
        self.depth_bins = depth_bins
        if not bkpt_bins:
            self.bkpt_bins = [depth_bins]
        else:
            self.bkpt_bins = bkpt_bins
        self.mapQT = mapq_t
        self.hcT = hardclip_t

        # initialise data structures
        self.depths = np.zeros((self.depth_bins.num, len(SamStats.depth_cols)), dtype=np.intc)
        self.aln_stats = []
        self.fwd_inserts = []
        self.rvs_inserts = []
        for bin in self.bkpt_bins:
            self.aln_stats.append(np.zeros((bin.num, len(SamStats.aln_stats_cols)), dtype=np.intc))
            self.fwd_inserts.append(np.empty(bin.num, dtype=list))
            self.rvs_inserts.append(np.empty(bin.num, dtype=list))
            for j in range(0, bin.num):
                self.fwd_inserts[-1][j] = []
                self.rvs_inserts[-1][j] = []

    def add_to_depth(self, coverage, cols):
        # add start and end covered bins, partial coverage
        for j in cols:
            self.depths[coverage[0][0]][j] += coverage[0][1]
            self.depths[coverage[1][0]][j] += coverage[1][1]

        # add in all bins that have full coverage
        for i in range(coverage[0][0]+1, coverage[1][0]):
            for j in cols:
                self.depths[i][j] += self.depth_bins.size

    def add_to_aln_stats(self, idx, coverage, cols):
        for i in range(coverage[0][0], coverage[1][0]+1):
            for j in cols:
                self.aln_stats[idx][i][j] += 1

    def process(self, sam_entry):
        cov = self.depth_bins.get_bin_coverage(sam_entry.pos, sam_entry.get_right_pos())
        if cov is None:
            return None
        depth_cols = [SamStats.TOTAL]

        if sam_entry.mapQ <= self.mapQT:
            if sam_entry.mapQ == 0:
                depth_cols.append(SamStats.MAPQ0)
            else:
                depth_cols.append(SamStats.MAPQLTT)

        self.add_to_depth(cov, depth_cols)

        for idx, bins in enumerate(self.bkpt_bins):
            cov = bins.get_bin_coverage(sam_entry.pos, sam_entry.get_right_pos())
            if cov is None:
                continue
            stats_cols = [SamStats.READS]

            if (sam_entry.flag & SamEntry.secondary):
                stats_cols.append(SamStats.SECONDARY)

            if (sam_entry.flag & SamEntry.supplementary):
                stats_cols.append(SamStats.SUPPLEMENTARY)

            if max(sam_entry.get_hard_clips()) >= self.hcT:
                stats_cols.append(SamStats.HARDCLIPPED)

            if sam_entry.has_unmapped_mate():
                stats_cols.append(SamStats.ORPHANED)
            else:
                if sam_entry.mate_same_strand():
                    stats_cols.append(SamStats.SAMESTRAND)
                elif sam_entry.is_inverted():
                    stats_cols.append(SamStats.INVERTED)
                else:
                    # correctly oriented pair reads
                    # filter low mapQ reads as these give spurious mapping distances
                    if not sam_entry.tlen == 0 and sam_entry.mapQ > self.mapQT:
                        # pairs are mapped correctly, so add to insert sizes
                        if sam_entry.is_rvs():
                            ins_cov = bins.get_bin_coverage(sam_entry.get_right_pos(), sam_entry.get_right_pos())
                            if ins_cov is None:
                                continue
                            self.rvs_inserts[idx][ins_cov[1][0]].append(-1 * sam_entry.tlen)
                        else:
                            ins_cov = bins.get_bin_coverage(sam_entry.pos, sam_entry.pos)
                            if ins_cov is None:
                                continue
                            self.fwd_inserts[idx][ins_cov[0][0]].append(sam_entry.tlen)
            self.add_to_aln_stats(idx, cov, stats_cols)

    def convert_depths(self):
        for i in range(0, self.depth_bins.num):
            for j in range(0, len(SamStats.depth_cols)):
                self.depths[i][j] /=  self.depth_bins.size

    # Print the collected stats
    def print_stats(self, dir):
        depth_file = file(os.path.join(dir, 'depths.tsv'), 'w')
        aln_stats_file = file(os.path.join(dir, 'aln_stats.tsv'), 'w')
        fwd_ins_file = file(os.path.join(dir, 'fwd_ins.tsv'), 'w')
        rvs_ins_file = file(os.path.join(dir, 'rvs_ins.tsv'), 'w')

        self.convert_depths()
        # print depths
        depth_file.write('bin\t' + '\t'.join(SamStats.depth_cols) + '\n')
        for i, row in enumerate(self.depths):
            depth_file.write(str(self.depth_bins.start + i*self.depth_bins.size) + '\t')
            for j in range(0, len(row)-1):
                depth_file.write(str(row[j]) + '\t')
            depth_file.write(str(row[-1]) + '\n')
        depth_file.close()

        # print alignmet stats and insert sizes
        aln_stats_file.write('bin\t' + '\t'.join(SamStats.aln_stats_cols) + '\n')
        for i, bin in enumerate(self.bkpt_bins):
            for j, row in enumerate(self.aln_stats[i]):
                # aln_stats
                aln_stats_file.write(str(bin.start + j*bin.size) + '\t')
                for k in range(0,len(row)-1):
                    aln_stats_file.write(str(row[k]) + '\t')
                aln_stats_file.write(str(row[k]) + '\n')
                # fwd inserts
                fwd_ins_file.write(str(bin.start + j * bin.size) + '\t')
                if self.fwd_inserts[i][j]:
                    for k in range(0, len(self.fwd_inserts[i][j])-1):
                         fwd_ins_file.write(str(self.fwd_inserts[i][j][k]) + ',')
                    fwd_ins_file.write(str(self.fwd_inserts[i][j][-1]))
                fwd_ins_file.write('\n')
                # rvs inserts
                rvs_ins_file.write(str(bin.start + j * bin.size) + '\t')
                if self.rvs_inserts[i][j]:
                    for k in range(0, len(self.rvs_inserts[i][j])-1):
                        rvs_ins_file.write(str(self.rvs_inserts[i][j][k]) + ',')
                    rvs_ins_file.write(str(self.rvs_inserts[i][j][-1]))
                rvs_ins_file.write('\n')
        fwd_ins_file.close()
        rvs_ins_file.close()
        aln_stats_file.close()

    # returns a list of sam_stats corresponding to the list of bams given for this position
    @staticmethod
    def get_sam_stats(chrom, start, end, sams, num_bins=100, breakpoints=None):
        depth_bins = Bins(start, end, ideal_num_bins=num_bins)
        region = chrom + ':' + str(depth_bins.start) + '-' + str(depth_bins.end)
        bkpt_bins = []
        # breakpoints is an iterable of regions (start, end)
        if breakpoints is not None:
            for bkpt in breakpoints:
                bkpt_bins.append(Bins(bkpt[0],bkpt[1], ideal_num_bins=num_bins/2))
        sam_stats = []
        for s in sams:
            sam_stats.append(SamStats(depth_bins, bkpt_bins=bkpt_bins))
            p = SAMtools.view(s, region)
            line = p.stdout.readline()
            while line:
                sam_stats[-1].process(SamEntry(line.split()))
                line = p.stdout.readline()
        return sam_stats


class SAMtools:
    @staticmethod
    def check_installation():
        cmd = ['samtools', '--version-only']
        try:
            subprocess.check_output(cmd)
        except OSError:
            print 'Error: could not run samtools. Are you sure it is installed?'
            exit(1)

    @staticmethod
    def view(sam, region, include_flag=None, exclude_flag=
            (SamEntry.duplicate + SamEntry.fails_QC + SamEntry.read_unmapped), samtools='samtools'):
        cmd = []
        cmd.append(samtools)
        cmd.append('view')
        if include_flag:
            cmd.append('-f')
            cmd.append(str(include_flag))
        if exclude_flag:
            cmd.append('-F')
            cmd.append(str(exclude_flag))
        cmd.append(sam)
        cmd.append(region)
        print ' '.join(cmd) + '\n'
        p = subprocess.Popen(cmd, bufsize=1024, stdout=PIPE)
        if p.poll():
            print "Error code %d from command:\n%s\n" % (' '.join(cmd) + '\n')
            exit(1)
        return p
