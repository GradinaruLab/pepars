from Aligner import Aligner
from sequencing import FASTQ
from workspace import Workspace as ws
import subprocess
import re
import os

import sys

class Bowtie_Aligner(Aligner):

    def __init__(self):
        Aligner.__init__(self)
        self.method = 'Bowtie'
        self.ids = []
        self.sequences = []

    def align_library(self, library, template, 
        allowed = ['G', 'T'],
        variants = [3,3,3,3], 
        sam_end = '_alignment.sam',
        cigar_end = '_bowtie_alignment_CIGAR_sequences.txt',
        reference_file = 'reference.fa', 
        reference_name = 'YichengReference'):
        """
        Aligns given library with the given template using Shashank's bottom-up
        method with Bowtie2
        """

        num_sequences = 0.0

        # Create references for bowtie to work with
        reference_name = ws.get_raw_data_path(reference_name)
        stringer = 'bowtie2-build ' + ws.get_raw_data_path(reference_file) + \
            ' '+ reference_name
        subprocess.call([stringer], shell = True, executable = '/bin/bash')

        reference = ws.get_raw_data_path(reference_file)

        for fastq_file_name in library.fastq_files:
            file_id = fastq_file_name.split('_')[0]
            fastq_file_name = ws.get_raw_data_path(fastq_file_name)
            sam_path = ws.get_raw_data_path(file_id + sam_end)

            # Create .sam file if not exists
            if not os.path.exists(sam_path):
                print sam_path
                cmd = ' '.join(['bowtie2 --local -x', reference_name,
                 '-U', fastq_file_name, '-S', ws.get_raw_data_path(file_id + sam_end)])
                subprocess.call([cmd], shell = True, executable = '/bin/bash')
            
            # From .sam file, create .txt file
            cmd = ' '.join(['cut -f1,6,10 ', sam_path,
                '> ' + ws.get_raw_data_path(file_id + cigar_end)])
            subprocess.call([cmd], shell=True, executable='/bin/bash')
            file = ws.get_raw_data_path(file_id+cigar_end)

            # Process the .txt file
            file = open(file, 'r').readlines()[3:]
            for line in file:
                line = line.split('\t')
                id = line[0]
                cigar = line[1]
                sequence = line[2]
                self.parse(id, cigar, sequence, allowed, variants)

            num_sequences += len(file)
        
        statistics = {}

        statistics['Number of Sequences'] = num_sequences
        statistics['Alignment rate'] = float(len(self.sequences))/float(num_sequences)

        return self.sequences, self.ids, statistics

    def parse(self, id, cigar, seq, allowed, variants):
        num_insertions = len(variants)
        index = 0
        chars = re.findall('[^\d]+', cigar)
        nums = re.findall('\d+', cigar)
        if nums.count('3') < num_insertions or chars.count('I') != num_insertions or \
            len(chars) != len(nums):
            return

        row = []
        for char, num in zip(chars, nums):
            num = int(num)
            if char=='I' and num == 3:
                codon = seq[index : index + num]
                if codon[-1] not in allowed:
                    return
                row.append(codon)
            index += num

        if len(row) != num_insertions:
            return
        # variants.append(row)
        self.sequences.append(seq)
        self.ids.append(id)