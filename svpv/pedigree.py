from __future__ import print_function

# read a pedigree and store info in a dictionary
class Pedigree:
    def __init__(self, ped, samples):
        # dict of list of samples in family
        self.samples_by_family = {}
        # dict of family by sample
        self.families_by_sample = {}
        # all info from ped file
        self.samples = {}
        for line in ped:
            if line[0] == '#' or line == '\n':
                continue
            try:
                fam, iid, fid, mid = line.split('\t')[0:4]
            except ValueError:
                print('Line in pedigree incorrectly formatted:\n"{}"\n'.format(line))
                continue
            if iid not in samples:
                continue
            if fam in self.samples_by_family:
                self.samples_by_family[fam].append(iid)
            else:
                self.samples_by_family[fam] = [iid]
            self.families_by_sample[iid] = fam
            self.samples[iid] = (fam, iid, fid, mid)

        for s in samples:
            if s not in self.samples:
                print('Warning: sample {} not in pedigree\n'.format(s))
