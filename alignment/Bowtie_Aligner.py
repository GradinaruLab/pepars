from Aligner import Aligner
from sequencing import FASTQ
from workspace import Workspace as ws
import subprocess
import os

class Bowtie_Aligner(Aligner):

    def __init__(self):
        Aligner.__init__(self)
        self.method = 'Bowtie'
        self.ids = []
        self.sequences = []

    def align_library(self, library, template, 
        num_sequences = 0.0,
        allowed = ['G', 'T'], variants = [3,3,3,3], sam_end = '_alignment.sam',
        cigar_end = '_bowtie_alignment_CIGAR_sequences.txt',
        reference_file = 'reference.fa', reference_name = 'YichengReference'):

        subprocess.call(['bowtie2-build', ws.get_raw_data_path(reference_file), reference_name])

        reference = ws.get_raw_data_path(reference_file) # put reference.fa file in raw_data
        for file_name in library.fastq_files:
            file_id = file_name.split('_')[0]
            file_name = ws.get_raw_data_path(file_name)
            sam_path = ws.get_raw_data_path(file_id + sam_end) # put .sam file in raw_data

            # Create .sam file if not exists
            if not os.path.exists(sam_path):
                print sam_path
                # subprocess.call(['bowtie2 --local -x', reference_name,
                #  '-U', file_name, '-S', file_id + sam_end])
            
            # From .sam file, create .txt file
            subprocess.call(['cut -f1,6,10 ', sam_path,
                '> ' + file_id + cigar_end])
            file = ws.get_raw_data_path(file_id+cigar_end)
            # Process the .txt file
            file = open(file, 'r').readlines()[3:]
            for line in file:
                line = line.split('\t')
                id = line[0]
                cigar = line[1]
                sequence = line[2]
                self.parse(id, cigar, sequence)

            num_sequences += len(file)
        
        statistics = {}

        statistics['Number of Sequences'] = num_sequences
        statistics['Alignment rate'] = float(len(self.sequences))/float(num_sequences)

        return self.sequences, self.ids, statistics

    def parse(id, cigar, seq):
        index = 0
        chars = re.findall('[^\d]+', cigar)
        nums = re.findall('\d+', cigar)
        if nums.count('3') < number_of_i or chars.count('I') != number_of_i or \
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

        if len(row) != number_of_i:
            return
        variants.append(row)
        self.sequences.append(seq)
        self.ids.append(id)