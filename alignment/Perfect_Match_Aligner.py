from .Aligner import Aligner
from sequencing import FASTQ
from workspace import Workspace as ws
from utils import DNA

class Perfect_Match_Aligner(Aligner):

    def __init__(self):
        Aligner.__init__(self)
        self.method = self.__class__.__name__

    def align_library(self, library, template, \
        variant_sequence_quality_threshold = 0, \
        mismatch_quality_threshold = 0,
        output_frequency = 1e5):

        self.template_sequence = template.sequence
        self.variant_sequence_quality_threshold = \
            int(variant_sequence_quality_threshold)
        self.mismatch_quality_threshold = int(mismatch_quality_threshold)
        template_mismatches = 0
        variant_quality_mismatches = 0
        variant_nucleotide_mismatches = 0
        size_mismatches = 0
        invalid_nucleotides = 0

        print('Aligning library \'' + library.name + '\'')

        template_length = len(self.template_sequence)

        extracted_sequences = []
        uuids = []
        error_probabilities = []
        sequence_counts = {}

        num_sequences = 0

        for fastq_file_name in library.fastq_files:

            fastq_file = ws.get_fastq_file(fastq_file_name).read().splitlines()

            line_count = 0

            for line in fastq_file:

                if line_count % 4 == 1:
                    if int(line_count / 4) % output_frequency == 0:
                        self.update_num_sequences_aligned(num_sequences)
                    sequence = line

                elif line_count % 4 == 3:

                    quality_string = line
                    
                    if (len(sequence)) == template_length:

                        extracted_sequence, uuid, error_type, error_probability = \
                            self.extract_sequence(sequence, quality_string)
                        
                        if error_type == 0:
                            extracted_sequences.append(extracted_sequence)
                            if extracted_sequence in sequence_counts:
                                sequence_counts[extracted_sequence] += 1
                            else:
                                sequence_counts[extracted_sequence] = 1
                            uuids.append(uuid)
                            error_probabilities.append(error_probability)
                        elif error_type == 1:
                            template_mismatches += 1
                        elif error_type == 2:
                            variant_quality_mismatches += 1
                        elif error_type == 3:
                            variant_nucleotide_mismatches += 1
                        elif error_type == 4:
                            invalid_nucleotides += 1
                    else:
                        size_mismatches +=1

                    num_sequences += 1

                line_count += 1

            ws.close_fastq_file(fastq_file_name)

        mismatched_sequences = template_mismatches + variant_quality_mismatches\
            + variant_nucleotide_mismatches + size_mismatches

        num_single_counts = 0
        for sequence, count in sequence_counts.items():
            if count == 1:
                num_single_counts += 1
        
        statistics = {}
        
        statistics["Number of Sequences"] = num_sequences
        statistics["Alignment rate"] = float(1.0-float(mismatched_sequences)/float(num_sequences))
        statistics["Template Mismatch Failure Rate"] = float(template_mismatches)/float(num_sequences)
        statistics["Variant Quality Failure Rate"] = float(variant_quality_mismatches)/float(num_sequences)
        statistics["Variant Nucleotide Mismatch Rate"] = float(variant_nucleotide_mismatches)/float(num_sequences)
        statistics["Template Size Mismatch Rate"] = float(size_mismatches)/float(num_sequences)
        statistics["Invalid nucleotide Rate"] = float(invalid_nucleotides)/float(num_sequences)
        statistics["Expected Number of Misreads"] = sum(error_probabilities)
        statistics["Number Single Count Sequences"] = num_single_counts

        return extracted_sequences, uuids, statistics

    def extract_sequence(self, fastq_line, quality_string):

        extracted_sequence = []
        uuid = []
        num_errors = 0
        probability_of_perfect_read = 1.0

        template_length = len(self.template_sequence)

        quality_score = FASTQ.convert_quality_string_to_quality_score(quality_string)
        #print("Quality score is: ")
        #print(quality_score)

        for template_idx in range(template_length):

            # Case 1: the template is an X, so we ignore this
            if self.template_sequence[template_idx] == 'X':
                # Do nothing
                continue
            # Case 2: The template is looking for a specific nudleotide, so we compare the 
            # template if and only if the quality score is above our threshold
            elif self.template_sequence[template_idx] in DNA.get_nucleotides():

                # If the quality score is below our threshold, we ignore this comparison
                if quality_score[template_idx] < self.mismatch_quality_threshold:
                    continue

                # If this isn't a match, it's an error
                if self.template_sequence[template_idx].upper() != fastq_line[template_idx].upper():
                    return False, False, 1, 1
            # Case 3: If this is not a valid nucleotide, it's an error
            elif fastq_line[template_idx] not in DNA.get_nucleotides():
                return False, False, 4, 1
            # Case 4: The template is an 'I', so this is part of the UUID
            elif self.template_sequence[template_idx] == 'I':
                uuid.append(fastq_line[template_idx])

            # Case 5: the template is an N, so this is part of our variant sequence
            else:
                if self.template_sequence[template_idx] == 'N':

                    if quality_score[template_idx] < int(self.variant_sequence_quality_threshold):
                        return False, False, 2, 1

                # Case 6: the template is a K, so the variant should be a T or G
                elif self.template_sequence[template_idx] == 'K':

                    if quality_score[template_idx] < self.variant_sequence_quality_threshold:
                        return False, False, 2, 1

                    if fastq_line[template_idx] != 'G' and fastq_line[template_idx] != 'T':
                        return False, False, 3, 1

                probability_of_error = FASTQ.convert_quality_score_to_probability_of_error(quality_score[template_idx])

                extracted_sequence.append(fastq_line[template_idx])
                probability_of_perfect_read *= (1.0 - probability_of_error)

                #print("Probability of perfect read: %.6f" % probability_of_perfect_read)

        return ''.join(extracted_sequence), ''.join(uuid), 0, 1.0 - probability_of_perfect_read