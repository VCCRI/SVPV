# # -*- coding: utf-8 -*-
# """
# author: Jacob Munro, Victor Chang Cardiac Research Institute
# """

import subprocess
from subprocess import PIPE
import copy


class VCFManager:
    vcf_count = 1

    def __init__(self, vcf_file, name='VCF ' + str(vcf_count), db_mode=False):
        VCFManager.vcf_count += 1
        self.name = name
        self.samples = BCFtools.get_samples(vcf_file)
        # count of svs in vcf
        self.count = 0
        #dict of Structural variants
        #chromosoms as keys
        #SVs['chr1'] = [list of SVs on chromosome 1]
        self.SVs = {}
        self.set_svs(vcf_file, db_mode)

    # for all svs, remove those that have not been called in the list of samples
    def remove_absent_svs(self, samples):
        idxs = []
        for s in samples:
            if s in self.samples:
                idxs.append(self.samples.index(s))

        for chrom in self.SVs:
            delete = []
            for i, sv in enumerate(self.SVs[chrom]):
                present = False
                for idx in idxs:
                    if '1' in sv.GTs[idx]:
                        present = True
                        break
                if not present:
                    delete.append(i)
            for i in sorted(delete, reverse=True):
                del self.SVs[chrom][i]
                self.count -= 1

    # read in all SV sites
    def set_svs(self, vcf, db_mode):
        p = BCFtools.get_SV_sites(vcf, db_mode)
        line = p.stdout.readline()
        while line:
            sv = SV.attempt_sv_parse(line.split(), db_mode)
            if sv is not None:
                self.count += 1
                if sv.chrom in self.SVs:
                    self.SVs[sv.chrom].append(sv)
                else:
                    self.SVs[sv.chrom] = [sv]
            line = p.stdout.readline()

    # return all SV calls that overlap with given range
    def get_svs_in_range(self, chrom, start, end, sample=None):
        ret = []
        if chrom in self.SVs:
            for sv in self.SVs[chrom]:
                if sv.start > end:
                    continue
                if sv.end < start:
                    continue
                else:
                    ret.append(sv)

        if sample is not None:
            delete = []
            for i, sv in enumerate(ret):
                if '1' not in sv.GTs[self.samples.index(sample)]:
                    delete.append(i)
            for i in sorted(delete, reverse=True):
                del ret[i]

        return ret

    def get_sample_index(self, sample):
        if sample in self.samples:
            return self.samples.index(sample)
    
    def filter_svs(self, filter_par, as_list=False):
        matching = {}
        if filter_par.chrom:
            if filter_par.chrom in self.SVs:
                matching[filter_par.chrom] = copy.copy(self.SVs[filter_par.chrom])
            else:
                # no SVs in chrom, return empty dict
                return matching
        else:
            matching = {}
            for chrom in self.SVs:
                matching[chrom] = copy.copy(self.SVs[chrom])

        # filter by svtype
        if filter_par.svtype:
            for chrom in matching:
                delete = []
                for i in range(0, len(matching[chrom])):
                    if matching[chrom][i].svtype != filter_par.svtype:
                        delete.append(i)
                delete.sort(reverse=True)
                for i in delete:
                    del matching[chrom][i]

        if filter_par.sample_GTs:
            # go through all SVs currently in matching and remove those without the correct genotype
            for chrom in matching:
                delete = []
                for i in range(0, len(matching[chrom])):
                    for sample in filter_par.sample_GTs:
                        if matching[chrom][i].GTs[self.get_sample_index(sample)] not in filter_par.sample_GTs[sample] \
                                and '*' not in filter_par.sample_GTs[sample]:
                            delete.append(i)
                            break
                delete.sort(reverse=True)
                for i in delete:
                    del matching[chrom][i]
        # filter by maf
        if filter_par.AF_thresh:
            for chrom in matching:
                delete = []
                for i, sv in enumerate(matching[chrom]):
                    if filter_par.AF_thresh_is_LT:
                        if not sv.get_AF() < filter_par.AF_thresh:
                            delete.append(i)
                    else:
                        if not sv.get_AF() > filter_par.AF_thresh:
                            delete.append(i)
                delete.sort(reverse=True)
                for i in delete:
                    del matching[chrom][i]

        # filter by ref_genes/intersection with specific gene
        if filter_par.RG_intersection or filter_par.gene_list_intersection or filter_par.exonic:
            for chrom in matching:
                delete = []
                for i, sv in enumerate(matching[chrom]):
                    intersecting = filter_par.ref_genes.get_entries_in_range(chrom, sv.start, sv.end)
                    if not intersecting:
                        delete.append(i)
                    elif filter_par.gene_list_intersection:
                        in_gene_list = False
                        for gene in intersecting:
                            if gene.name2.upper() in filter_par.gene_list:
                                if filter_par.exonic and not gene.intersects_exon(sv.start, sv.end):
                                    continue
                                in_gene_list = True
                                break
                        if not in_gene_list:
                            delete.append(i)
                    elif filter_par.exonic:
                        exonic = False
                        for gene in intersecting:
                            if gene.intersects_exon(sv.start, sv.end):
                                exonic = True
                                break
                        if not exonic:
                            delete.append(i)

                delete.sort(reverse=True)
                for i in delete:
                    del matching[chrom][i]

        # filter by SV length
        if (not filter_par.min_len == None) or (not filter_par.max_len == None):
            for chrom in matching:
                delete = []
                for i, sv in enumerate(matching[chrom]):
                    if not filter_par.min_len == None:
                        if (sv.end - sv.start + 1) < filter_par.min_len:
                            delete.append(i)
                            continue
                    if not filter_par.max_len == None:
                        if (sv.end - sv.start + 1) > filter_par.max_len:
                            delete.append(i)
                            continue
                delete.sort(reverse=True)
                for i in delete:
                    del matching[chrom][i]

        # return the filtered set of SVs
        if as_list:
            listed = []
            for chrom in sorted(matching.keys()):
                listed.extend(matching[chrom])
            return listed
        return matching


class SV:
    valid_SVs = ['DEL', 'DUP', 'CNV', 'INV']

    @staticmethod
    def attempt_sv_parse(fields, db_mode):
        if not len(fields) > 3 or not fields[3] in SV.valid_SVs:
            return None
        else:
            return SV(fields, db_mode)

    def __init__(self, fields, db_mode):
        self.chrom = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        self.svtype = fields[3]
        if db_mode:
            self.GTs = None
            ''' TODO '''
            ''' need to do this properly for multi-allelic positions - currently just showing allele B freq '''
            self.AF = float(fields[4].split(',')[0])
        else:
            self.GTs = fields[4:]
            self.AF = self.get_AF()

    def get_AF(self):
        if len(self.GTs):
            n = 0
            for gt in self.GTs:
                if '1' in gt:
                    if '1/1' in gt:
                        n += 2
                    else:
                        n += 1
            return (n / float(2 * len(self.GTs)))

    # helper method for print_SVs
    def to_string(self, sample_index=None):
        if not sample_index == None:
            return '\t'.join([self.chrom, str(self.start), str(self.end), self.svtype, self.GTs[sample_index]]) + '\n'
        else:
            return '\t'.join([self.chrom, str(self.start), str(self.end), self.svtype, str(self.AF)]) + '\n'

    def string_tuple(self):
        return (self.svtype, self.chrom, str(self.start), str(self.end - self.start +1),
                ('{0:.2f}'.format(self.get_AF())))

    # print SVs, either with genotype per sample, or MAF for whole BATCH annotation
    @staticmethod
    def print_SVs(SVs, out, vcf_name, sample_index=None):
        for sv in SVs:
            out.write(vcf_name + '\t' + sv.to_string(sample_index=sample_index))

    @staticmethod
    def print_SVs_header(out, sample_index=None):
        if not sample_index == None:
            out.write('\t'.join(['vcf', 'chrom', 'start', 'end', 'svtype', 'gt']) + '\n')
        else:
            out.write('\t'.join(['vcf', 'chrom', 'start', 'end', 'svtype', 'MAF']) + '\n')


class BCFtools:
    @staticmethod
    def check_installation():
        cmd = ['bcftools']
        try:
            subprocess.Popen(cmd)
        except OSError:
            print 'Error: could not run bcftools. Are you sure it is installed?'
            exit(1)

    # return a pipe to the set of sv sites
    @staticmethod
    def get_SV_sites(vcf, db_mode=False):
        cmd = []
        cmd.append("bcftools")
        cmd.append("query")
        cmd.append("-f")
        if db_mode:
            cmd.append("%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%INFO/AF\\n")
        else:
            cmd.append("%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE[\\t%GT]\\n")
        cmd.append(vcf)
        print ' '.join(cmd) + '\n'
        p = subprocess.Popen(cmd, bufsize=1024, stdout=PIPE)
        if p.poll():
            print "Error code %d from command:\n%s\n" % (' '.join(cmd) + '\n')
            exit(1)
        return p

    #return a list of samples
    @staticmethod
    def get_samples(vcf):
        cmd = []
        cmd.append("bcftools")
        cmd.append("query")
        cmd.append("-l")
        cmd.append(vcf)
        print ' '.join(cmd) + '\n'
        return subprocess.check_output(cmd).split()
