import subprocess
import re
import os

from .Aligner import Aligner
from ..utils import DNA


class Bowtie_Aligner(Aligner):

    def __init__(self):
        Aligner.__init__(self)
        self.method = 'Bowtie'
        self.ids = []
        self.sequences = []

    def _align(
            self, FASTQ_file_sets, template,
            reverse_complement_template,
            working_directory=os.getcwd(),
            is_local = False,
            output_frequency = 1e5,
            approach = None,
            allow_insertions_deletions = False,
            quality_threshold = 0.0):
        """
        Aligns given library with the given template using Shashank's bottom-up
        method with Bowtie2
        """

        ##
        reference_name = '_'.join(['bowtie', approach])
        other_name = '_'.join(['bowtie',approach])
        
        ##
        reference_file_name = reference_name + '.fa'

        if approach == 'elimination':
            parse_function = self.parse_elimination
        elif approach == 'mismatch':
            parse_function = self.parse_nnt

        variants = self.get_variants(template)

        num_sequences = 0.0

        reference_file_name = \
            os.path.join(working_directory, reference_file_name)
        # Create reference file if not exists
        try:
            reference_file = open(reference_file_name, 'r')
        except:
            reference_file = open(reference_file_name, 'w')

            if approach == 'elimination':
                template.replace('NNK','')
            elif approach == 'mismatch':
                template.replace('K','T')
            
            reference_file.write('>' + reference_name + os.linesep)
            reference_file.write(template)
        
        reference_file.close()

        # Create references for bowtie to work with
        stringer = 'bowtie2-build ' + reference_file_name + \
            ' '+ reference_name

        subprocess.call([stringer], shell = True, executable = '/bin/bash')

        for FASTQ_file_set in FASTQ_file_sets:

            if (num_sequences / 4) % output_frequency == 0:
                self.update_num_sequences_aligned(num_sequences)
            
            # file_id = fastq_file_name.split('_')[0]

            fastq_file_path = os.path.join(working_directory,
                                           FASTQ_file_set.files[0].file_name)
            sam_path = os.path.join(working_directory, other_name + '.sam')

            # Create .sam file if not exists
            if not os.path.exists(sam_path):
                if not is_local:
                    cmd = ' '.join(['bowtie2 -x', reference_name,
                        '-U', fastq_file_path, '-S',
                        sam_path])
                else:
                    cmd = ' '.join(['bowtie2 --local -x', reference_name,
                        '-U', fastq_file_path, '-S',
                        sam_path])

                subprocess.call([cmd], shell = True, executable = '/bin/bash')
            
            # From .sam file, create .txt file
            cmd = ' '.join(['cut -f1,6,10,11 ', sam_path,
                '> ' + os.path.join(working_directory, other_name + '.txt')])
            subprocess.call([cmd], shell=True, executable='/bin/bash')
            file = os.path.join(working_directory, other_name + '.txt')

            # Process the .txt file
            file = open(file, 'r').readlines()[3:]
            for line in file:
                line = line.split('\t')
                id = line[0]
                cigar = line[1]
                sequence = line[2]
                quality = line[3]
                parse_function(id, cigar, sequence, variants)

            num_sequences += len(file)
        
        statistics = {}

        statistics['Number of Sequences'] = num_sequences
        statistics['Alignment rate'] = float(len(self.sequences))/\
                                       float(num_sequences)

        return self.sequences, self.ids, statistics


    def parse_nnt(self, id, cigar, seq, variants):
        pass

    def parse_elimination(self, id, cigar, seq, variants):
        num_insertions = len(variants)
        index = 0
        chars = re.findall('[^\d]+', cigar)
        nums = re.findall('\d+', cigar)
        
        if chars.count('I') != num_insertions or \
            len(chars) != len(nums):
            return

        row = []
        for i in range(len(chars)):
            
            char = chars[i]
            num = int(nums[i])

            variant = variants[len(row)]
            if char=='I' and num == len(variant):
                codon = seq[index : index + num]
                if not self.match(codon, variant):
                    continue
                row.append(codon)
            index += num

        if len(row) != num_insertions:
            return
        # variants.append(row)

        self.sequences.append(seq)
        self.ids.append(id)

    def get_variants(self, sequence):
        # Add more IUPAC characters
        chars = ''.join(DNA.IUPAC_GRAMMAR_MAP.items())
        variants = re.findall('['+chars+']+',sequence)
        return variants

    def match(self, sequence_fraction, template_fraction):
        assert len(codon)==len(template_fraction)

        for i in range(len(template_fraction)):
            if sequence_fraction[i] in ('I', 'X'):
                return True
            if not sequence_fraction[i] in \
                   DNA.IUPAC_GRAMMAR_MAP[template_fraction[i]]:
                return False

        return True
