import numpy as np
import re
import time
from os import listdir
from os.path import isfile, join
import subprocess

codons = [3, 3, 3, 3]
num_insertions = len(codons)
allowed = ['G', 'T']

limit = 0
qlt_limit = 0.5

print subprocess.call('pwd')
names = {}
end = '_bowtie_alignment_CIGAR_sequences.txt'
my_path = '/home/peturhelgi/documents/code/virusfarm/Data/TestData/'
# files = [file for file in listdir(my_path) if isfile(join(my_path,file)) and end in file]
#files = ['17124_bowtie_alignment_CIGAR_sequences.txt']

##########################
# files=['17124_bowtie_alignment_CIGAR_sequences.txt']
# files=['example_17124_cigar.txt']
files=['quality_test.txt']


'''
 1456  bowtie2-build reference.fa YichengReference

 1459  bowtie2 --local -x YichengReference -U 17124_AGTTCC_L001_R1_001.fastq -S 17124_alignment0804.sam
 1462  cut -f1,6,10 17124_alignment0804.sam > 17124_cigar0804.txt
'''


def parse(id, cigar, seq):
    index = 0
    chars = re.findall('[^\d]+',cigar)
    nums = re.findall('\d+',cigar)
    index_of_var = 0
    indices = []
#####################################################

    # print chars
    # print nums

#####################################################


    row = []
    for i, (char, num) in enumerate(zip(chars, nums)):
        num = int(num)
        if len(row) >= len(codons):
            break
        if char=='I' and num == 3:
            codon = seq[index:index+num]
            index_of_var +=1
            if codon[-1] not in allowed:
                continue
                # pass
            indices.append(index_of_var)
            row.append(codon)
        if char != 'D':
            index += num

######################################################

    # print row, indices

######################################################

    if len(row) != num_insertions:
        return
    # if chars == ['S', 'M', 'I', 'M', 'I', 'M', 'I', 'M', 'I', 'M', 'I', 'M', 'S']:
    #     print seq
    variants.append(row)
    seqs.append(seq)
    ids.append(id)

count = 0
for file in files:
    print my_path,file
    ids = []
    variants = []
    seqs = []
    name = file
    start = time.clock()
    if limit != 0:
        file = open(join(my_path,file), 'r').readlines()[3:limit+3]
    else:
        print 'limit'
        file = open(join(my_path,file), 'r').readlines()[3:]
        limit = len(file)
        
    for line in file:
        line = line.split('\t')
        id = line[0]
        cigar = line[1]
        seq = line[2]
        qlt = line[3]
        if float(len(re.findall('\w', qlt)))/float(len(qlt)) < qlt_limit:
            count += 1
            continue
        if cigar!='*':
            parse(id, cigar, seq)
    end = time.clock()
    print name, 'took', end-start,'seconds'

    if not (len(ids) == len(seqs) == len(variants)):
        print 'Not consistent'
    print 'File had', len(variants), '/', len(file)-count,'matches'
    acc = (0.0+len(variants))/(0.0+len(file)-count)
    print 'Accuracy:', acc,'\n'
    print 'Number of seqs deleted due to low quality:',count,'out of' ,limit, '. That is',str(float(count)/float(limit))+'%'
    names[name]={'Total':len(file), 'Matched':len(ids), 'Accuracy':acc}

for name in sorted(names):
    print name
