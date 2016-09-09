# # -*- coding: utf-8 -*-
# """
# @author: j.munro@victorchang.edu.au
#
# """


class RefgeneManager:
    def __init__(self, ref_genes, keep_all=False):
        # dict by chrom (as with SVs in VCF_Manager)
        self.entries = {}
        for line in file(ref_genes):
            if line[0] == '#':
                continue
            fields = line.split()
            if not len(fields) == 16:
                print "Incorrect Number of Fields for RefGene Schema"
                return None
            e = RefGeneEntry(fields, keep_all=keep_all)
            if e.chrom in self.entries:
                self.entries[e.chrom].append(e)
            else:
                self.entries[e.chrom] = [e]

    def get_entries_in_range(self, chrom, start, end):
        ret = []
        if chrom in self.entries:
            for e in self.entries[chrom]:
                if e.txStart > end:
                    break
                if e.txEnd < start:
                    continue
                else:
                    ret.append(e)
        return ret


class RefGeneEntry:
    header = '\t'.join(['chrom', 'strand', 'txStart', 'txEnd', 'exonStarts', 'exonEnds', 'name2']) + '\n'

    def __init__(self, fields, keep_all=False):
        self.chrom = fields[RGF.chrom]
        self.strand = fields[RGF.strand]
        self.txStart = int(fields[RGF.txStart])
        self.txEnd = int(fields[RGF.txEnd])
        self.exonStarts = fields[RGF.exonStarts]
        self.exonEnds = fields[RGF.exonEnds]
        self.name2 = fields[RGF.name2]
        if keep_all:
            self.fields = fields

    def intersects_exon(self, other_start, other_end):
        starts = self.exonStarts.split(',')
        ends = self.exonEnds.split(',')
        for i in range(0, len(starts)):
            try:
                start = int(starts[i])
                end = int(ends[i])
            except ValueError:
                continue
            if other_start <= start <= other_end:
                return True
            elif end >= other_start > start:
                return True
        return False

    def to_string(self):
        return '\t'.join([self.chrom, self.strand, str(self.txStart), str(self.txEnd), self.exonStarts, self.exonEnds, self.name2]) + '\n'

    @staticmethod
    def print_entries(entries, out):
        out.write(RefGeneEntry.header)
        for e in entries:
            out.write(e.to_string())
        out.close()


# refGene Fields
class RGF:
    bin = 0
    name = 1
    chrom = 2
    strand = 3
    txStart = 4
    txEnd = 5
    cdsStart = 6
    cdsEnd = 7
    exonCount = 8
    exonStarts = 9
    exonEnds = 10
    score = 11
    name2 = 12
    cdsStartStat = 13
    cdsEndStat = 14
    exonFrames = 15
