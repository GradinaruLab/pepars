import numpy as np
import re
import time
from os import listdir
from os.path import isfile, join
codons = [3, 3, 3, 3]
num_insertions = len(codons)
# allowed = ['G', 'T', 'A', 'C']
allowed = ['G', 'T', 'A']
allowed.append('N')
'''
 1456  bowtie2-build reference.fa YichengReference

 1459  bowtie2 --local -x YichengReference -U 17115_ATCACG_L001_R1_001.fastq -S 17115_alignment.sam
 1462  cut -f1,6,10 17115_alignment.sam > 17115_bowtie_alignment_CIGAR_sequences.txt
'''
end = '_bowtie_alignment_CIGAR_sequences.txt'
my_path = '../Data/TestData/'

def parse(id, cigar, seq):
    index = 0
    chars = re.findall('[^\d]+',cigar)
    nums = re.findall('\d+',cigar)

    row = []
    for i, (char, num) in enumerate(zip(chars, nums)):
        num = int(num)
        if len(row) >= len(codons):
            break
        if char=='I' and num == 3:
            codon = seq[index:index+num]
            if codon[-1] not in allowed:
                continue
            row.append(codon)
        if char != 'D':
            index += num

    if len(row) != num_insertions:
        return
    variants.append(row)
    seqs.append(seq)
    ids.append(id)
# files = [file for file in listdir(my_path) if isfile(join(my_path,file)) and end in file]
files = ['17124_bowtie_alignment_CIGAR_sequences.txt']
names = {}
for file in files:
    ids = []
    variants = []
    seqs = []
    name = file
    start = time.clock()
    file = open(join(my_path,file), 'r').readlines()[3:]

    for line in file:
        line = line.split('\t')
        id = line[0]
        cigar = line[1]
        seq = line[2]
        parse(id, cigar, seq)
    end = time.clock()
    print name, 'took', end-start,'seconds'

    if not (len(ids) == len(seqs) == len(variants)):
        print 'Not consistent'
    print 'File had', len(variants), '/', len(file),'matches'
    acc = (0.0+len(variants))/(0.0+len(file))
    print 'Accuracy:', acc,'\n'

    names[name]={'Total':len(file), 'Matched':len(ids), 'Accuracy':acc}

for name in sorted(names):
    print name
