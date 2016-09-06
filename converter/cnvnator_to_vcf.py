# # -*- coding: utf-8 -*-
# """
# @author: j.munro@victorchang.edu.au
#
# """

import sys
import copy
import re
import math
import os

def main():
    usage = "[Required args]\n" \
        "-samples\t/file/with/samples/and/call_file/paths.txt\n" \
        "-o\t\t/out/location/and/name.vcf\n" \
        "-thresh\tjaccard index threshold to use for clustering\n" \
        "-chroms\tcomma separated list of chromosomes to process. defaults to all.\n"

    sample_calls = None
    out_vcf = None
    chroms = []
    thresh = 0.7
    # parse arguments
    for i, arg in enumerate(sys.argv):
        if arg == '-samples':
            sample_calls = file(sys.argv[i+1], 'r')
        elif arg == '-o':
            out_vcf = file(sys.argv[i+1], 'w')
        elif arg == '-chroms':
            chroms = sys.argv[i+1].split(',')
        elif arg == '-thresh':
            thresh = float(sys.argv[i+1])
            if thresh < 0.5:
                print 'Please enter a threshold above 0.5'
                exit(1)

    if not (sample_calls and out_vcf):
        print usage
        exit(1)

    samples = parse_samples(sample_calls, chroms)
    all_calls, edges = get_all_calls(samples, thresh)
    clustered = cluster(all_calls, edges, thresh)
    print_vcf(samples, clustered, out_vcf)


def parse_samples(sample_calls, chroms):
    print 'parsing sample files'
    samples = []
    for line in sample_calls:
        try:
            sample_name, sample_file = line.split()
        except ValueError:
            print '%s skipped due to incorrect number of fields' % line
            continue
        if not os.path.isfile(sample_file):
            print 'file not found: %s' % sample_file
            continue
        samples.append(Sample(sample_name))
        samples[-1].parse_calls(sample_file, chroms)
    print('\tdone.\t\n')
    return samples


class Edges():
    def __init__(self):
        self.dict = {}

    def add(self, chrom, svtype, sv1, sv2):
        if chrom not in self.dict:
            self.dict[chrom] = {}
        if svtype not in self.dict[chrom]:
            self.dict[chrom][svtype] = {}

        if sv1 not in self.dict[chrom][svtype]:
            self.dict[chrom][svtype][sv1] = [sv2]
        elif sv2 not in self.dict[chrom][svtype][sv1]:
            self.dict[chrom][svtype][sv1].append(sv2)
        self.dict[chrom][svtype][sv1].sort(key=sv1.jaccard, reverse=True)

        if sv2 not in self.dict[chrom][svtype]:
            self.dict[chrom][svtype][sv2] = [sv1]
        elif sv1 not in self.dict[chrom][svtype][sv2]:
            self.dict[chrom][svtype][sv2].append(sv1)
        self.dict[chrom][svtype][sv2].sort(key=sv2.jaccard, reverse=True)

    def update(self, chrom, svtype, sv1, sv2, new):
        if new not in self.dict[chrom][svtype]:
            self.dict[chrom][svtype][new] = []

        for n in self.dict[chrom][svtype][sv1]:
            if n is sv1 or n is sv2:
                continue
            while sv1 in self.dict[chrom][svtype][n]:
                self.dict[chrom][svtype][n].remove(sv1)
            self.dict[chrom][svtype][n].append(new)
            self.dict[chrom][svtype][n].sort(key=n.jaccard, reverse=True)
            if n not in self.dict[chrom][svtype][new] and n is not new:
                self.dict[chrom][svtype][new].append(n)

        for n in self.dict[chrom][svtype][sv2]:
            if n is sv1 or n is sv2:
                continue
            while sv2 in self.dict[chrom][svtype][n]:
                self.dict[chrom][svtype][n].remove(sv2)
            self.dict[chrom][svtype][n].append(new)
            self.dict[chrom][svtype][n].sort(key=n.jaccard, reverse=True)
            if n not in self.dict[chrom][svtype][new] and n is not new:
                self.dict[chrom][svtype][new].append(n)

        while sv1 in self.dict[chrom][svtype][new]:
            self.dict[chrom][svtype][new].remove(sv1)
        while sv2 in self.dict[chrom][svtype][new]:
            self.dict[chrom][svtype][new].remove(sv2)
        while new in self.dict[chrom][svtype][new]:
            self.dict[chrom][svtype][new].remove(new)
        self.dict[chrom][svtype][new].sort(key=new.jaccard, reverse=True)


def get_all_calls(samples, thresh):
    print 'compiling all calls:'
    all_calls = {}
    for sample in samples:
        for chrom in sample.calls:
            if chrom not in all_calls:
                all_calls[chrom] = {}
            for svtype in sample.calls[chrom]:
                if svtype not in all_calls[chrom]:
                    all_calls[chrom][svtype] = {}
                for bin in sample.calls[chrom][svtype]:
                    if bin not in all_calls[chrom][svtype]:
                        all_calls[chrom][svtype][bin] = []
                    all_calls[chrom][svtype][bin].extend(sample.calls[chrom][svtype][bin])
                    all_calls[chrom][svtype][bin].sort(key=lambda x: x.start)
    print('\tdone.\t\n')
    # connect neighbours
    # this done by passing through and connecting all those that are within +/- 1 size bin of bin of call
    # adding only those with a start poisition >= our start position
    print 'preprocessing:'
    edges = Edges()
    for chrom in all_calls:
        for svtype in all_calls[chrom]:
            print('\t%s, %s\t' % (chrom, svtype))
            for bin in all_calls[chrom][svtype]:
                # need to check bin-1 and bin+1 for neighbours as well
                if (bin - 1) in all_calls[chrom][svtype]:
                    lt_idx = 0
                else:
                    lt_idx = None
                if (bin + 1) in all_calls[chrom][svtype]:
                    gt_idx = 0
                else:
                    gt_idx = None
                for i, sv in enumerate(all_calls[chrom][svtype][bin]):

                    j = i+1
                    # add neighbours from this size bin
                    while j < len(all_calls[chrom][svtype][bin]) and sv.overlaps(all_calls[chrom][svtype][bin][j]):
                        if sv.jaccard(all_calls[chrom][svtype][bin][j]) > thresh:
                            edges.add(chrom, svtype, sv, all_calls[chrom][svtype][bin][j])
                        j += 1
                    # add neighbours from lesser size bin
                    if lt_idx is not None:
                        while lt_idx < len(all_calls[chrom][svtype][bin-1]):
                            if sv.start <= all_calls[chrom][svtype][bin-1][lt_idx].start:
                                break
                            lt_idx += 1

                        idx = lt_idx
                        while (idx < len(all_calls[chrom][svtype][bin-1]) and
                               sv.overlaps(all_calls[chrom][svtype][bin-1][idx])):
                            if sv.jaccard(all_calls[chrom][svtype][bin-1][idx]) > thresh:
                                edges.add(chrom, svtype, sv, all_calls[chrom][svtype][bin-1][idx])
                            idx += 1

                    # add neighbours from greater size bin
                    if gt_idx is not None:
                        while gt_idx < len(all_calls[chrom][svtype][bin+1]):
                            if sv.start <= all_calls[chrom][svtype][bin+1][gt_idx].start:
                                break
                            gt_idx += 1

                        idx = gt_idx
                        while (idx < len(all_calls[chrom][svtype][bin+1]) and
                                sv.overlaps(all_calls[chrom][svtype][bin+1][idx])):
                            if sv.jaccard(all_calls[chrom][svtype][bin+1][idx]) > thresh:
                                edges.add(chrom, svtype, sv, all_calls[chrom][svtype][bin+1][idx])
                            idx += 1

    for chrom in all_calls:
        for svtype in all_calls[chrom]:
            merged_bins = []
            for bin in all_calls[chrom][svtype]:
                merged_bins.extend(all_calls[chrom][svtype][bin])
            all_calls[chrom][svtype] = merged_bins
    print('\tdone.\t\n')

    return all_calls, edges


# cluster sorted lists of all svs
def cluster(all_calls, edges, thresh):
    print 'clustering:'
    for chrom in all_calls:
        for svtype in all_calls[chrom]:
            print('\t%s, %s\t' % (chrom, svtype))
            progress = Cluster.get_progress()
            while True:
                print 'len clusters: %d, progress: %s' % (len(all_calls[chrom][svtype]), str(progress))
                for sv in all_calls[chrom][svtype]:
                    if sv not in edges.dict[chrom][svtype] or len(edges.dict[chrom][svtype][sv]) == 0:
                        continue
                    other = edges.dict[chrom][svtype][sv][0]
                    if other is sv:
                        print 'Error: self reference'
                        exit(1)
                    if sv.jaccard(other) > thresh and (sv is edges.dict[chrom][svtype][other][0] or sv.jaccard(other) >
                                                       other.jaccard(edges.dict[chrom][svtype][other][0]) - 1e6):
                        new = Cluster.cluster_pair(sv, other)
                        edges.update(chrom, svtype, sv, other, new)
                        all_calls[chrom][svtype].remove(other)
                        all_calls[chrom][svtype].remove(sv)
                        all_calls[chrom][svtype].append(new)
                if progress == Cluster.get_progress():
                    break
                else:
                    progress = Cluster.get_progress()
    return all_calls


# print out vcf entries
def print_vcf(samples, clustered, out_vcf, debug=False):
    if debug:
        count = 0
        for chrom in clustered:
            for svtype in clustered[chrom]:
                clustered[chrom][svtype].sort(key=lambda x: x.start)
                for i, sv in enumerate(clustered[chrom][svtype]):
                    j = i+1
                    if j >= len(clustered[chrom][svtype]):
                        break
                    other = clustered[chrom][svtype][j]
                    while sv.overlaps(other) and j < (len(clustered[chrom][svtype])-1):
                        if sv.jaccard(other) > 0.7:
                            count += 1
                            if count < 10:
                                print sv.jaccard(other)
                                print sv
                                print '%s,%d %d bp' % (sv.chrom, sv.start, (sv.end - sv.start))

                                print other
                                print '%s,%d %d bp\n' % (other.chrom, other.start, (other.end - other.start))
                        j += 1
                        other = clustered[chrom][svtype][j]
        print count

    else:
        sample_names = []
        default_gts = []
        for s in samples:
            sample_names.append(s.name)
            default_gts.append('0/0')

        out_vcf.write(Header + '\t'.join(sample_names) + '\n')

        for chrom in sorted(clustered.keys()):
            merged_types = []
            for svtype in clustered[chrom]:
                merged_types.extend(clustered[chrom][svtype])
            merged_types.sort(key=lambda x: x.start)
            for sv in merged_types:
                gts = copy.copy(default_gts)
                if isinstance(sv, Cluster):
                    for call in sv.members:
                        gts[sample_names.index(call.sample_name)] = call.genotype()
                else:
                    gts[sample_names.index(sv.sample_name)] = sv.genotype()
                out_vcf.write(sv.to_vcf() + '\t' + '\t'.join(gts) + '\n')


Header = '##fileformat=VCFv4.2\n' \
    '##ALT=<ID=DEL,Description="Deletion">\n' \
    '##ALT=<ID=DUP,Description="Duplication">\n' \
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n' \
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n' \
    '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\n' \
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' \
    '##contig=<ID=chr1,length=248956422>\n' \
    '##contig=<ID=chr10,length=133797422>\n' \
    '##contig=<ID=chr11,length=135086622>\n' \
    '##contig=<ID=chr12,length=133275309>\n' \
    '##contig=<ID=chr13,length=114364328>\n' \
    '##contig=<ID=chr14,length=107043718>\n' \
    '##contig=<ID=chr15,length=101991189>\n' \
    '##contig=<ID=chr16,length=90338345>\n' \
    '##contig=<ID=chr17,length=83257441>\n' \
    '##contig=<ID=chr18,length=80373285>\n' \
    '##contig=<ID=chr19,length=58617616>\n' \
    '##contig=<ID=chr2,length=242193529>\n' \
    '##contig=<ID=chr20,length=64444167>\n' \
    '##contig=<ID=chr21,length=46709983>\n' \
    '##contig=<ID=chr22,length=50818468>\n' \
    '##contig=<ID=chr3,length=198295559>\n' \
    '##contig=<ID=chr4,length=190214555>\n' \
    '##contig=<ID=chr5,length=181538259>\n' \
    '##contig=<ID=chr6,length=170805979>\n' \
    '##contig=<ID=chr7,length=159345973>\n' \
    '##contig=<ID=chr8,length=145138636>\n' \
    '##contig=<ID=chr9,length=138394717>\n' \
    '##contig=<ID=chrX,length=156040895>\n' \
    '##contig=<ID=chrY,length=57227415>\n' \
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'


class SV:
    def __init__(self, chrom, start, end, svtype, rd=None):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.svtype = svtype
        self.sample_name = None
        if rd:
            self.rd = float(rd)
        else:
            self.rd = None

    def jaccard(self, other):
        if not self.chrom == other.chrom:
            return float(0)
        if self.end > other.end:
            max_end = self.end
            min_end = other.end
        else:
            max_end = other.end
            min_end = self.end
        if self.start < other.start:
            min_start = self.start
            max_start = other.start
        else:
            min_start = other.start
            max_start = self.start
        return float(min_end - max_start + 1) / float(max_end - min_start + 1)

    def genotype(self):
        if self.rd is not None and self.svtype == 'DEL':
            if self.rd < 0.25:
                return '1/1'
            else:
                return '0/1'
        elif self.rd is not None and self.svtype == 'DUP':
            if self.rd > 1.8:
                return '1/1'
            else:
                return '0/1'
        else:
            return '0/1'

    def overlaps(self, other):
        return not (other.end < self.start or other.start > self.end)

    def to_vcf(self):
        return '\t'.join([self.chrom, str(self.start), '.', '.', '<'+self.svtype+'>', '.', '.', 'IMPRECISE;SVTYPE=%s;END=%d' % (self.svtype, self.end), 'GT'])


class Cluster(SV):
    merges = 0
    clusters = 0
    additions = 0

    @staticmethod
    def get_progress():
        return (Cluster.merges, Cluster.clusters, Cluster.additions)

    def __init__(self, r1, r2):
        SV.__init__(self, r1.chrom, (r1.start + r2.start) / 2, (r1.end + r2.end) / 2, r1.svtype)
        self.members = [r1, r2]

    def update(self):
        start_sum = 0
        end_sum = 0
        for mem in self.members:
            start_sum += mem.start
            end_sum += mem.end
        self.start = start_sum / len(self.members)
        self.end = end_sum / len(self.members)

    def add_to_cluster(self, region):
        self.members.append(region)
        self.update()

    def merge_clusters(self, other):
        self.members.extend(other.members)
        self.update()

    @staticmethod
    def cluster_pair(r1, r2):
        if isinstance(r1, Cluster):
            if isinstance(r2, Cluster):
                Cluster.merges += 1
                r1.merge_clusters(r2)
                return r1
            else:
                Cluster.additions += 1
                r1.add_to_cluster(r2)
                return r1
        else:
            if isinstance(r2, Cluster):
                Cluster.additions += 1
                r2.add_to_cluster(r1)
                return r2
            else:
                Cluster.clusters += 1
                return Cluster(r1, r2)


class Sample:
    def __init__(self, name):
        self.name = name
        # dict[chrom] of dict[svtype] of dict[bin] of lists of svs
        self.calls = {}

    def parse_calls(self, CNVnator_file, chroms, f1=1e-6, f2=1e-6, f3=1e-6, f4=1e-6 ):
        for line in file(CNVnator_file):
            (svtype, coor, length, read_depth, e1, e2, e3, e4) = line.split()[0:8]
            if float(e1) > f1 and float(e2) > f2 and float(e3) > f3 and float(e4) > f4:
                continue
            (chrom, start, end) =  re.split('[:-]', coor)
            if chroms and chrom not in chroms:
                continue
            sv = SV(chrom, start, end, re.sub('deletion', 'DEL', re.sub('duplication', 'DUP', svtype)), rd=float(read_depth))
            self.add_call(sv)

    def add_call(self, call):
        if call.chrom not in self.calls:
            self.calls[call.chrom] = {}
        if call.svtype not in self.calls[call.chrom]:
            self.calls[call.chrom][call.svtype] = {}
        bin = int(math.log((call.end - call.start), 2))
        if bin not in self.calls[call.chrom][call.svtype]:
            self.calls[call.chrom][call.svtype][bin] = []
        call.sample_name = self.name
        self.calls[call.chrom][call.svtype][bin].append(call)

if __name__ == "__main__":
    main()
