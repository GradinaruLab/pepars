from Aligner import Aligner
from sequencing import FASTQ
from workspace import Workspace as ws

class Perfect_Match_Aligner(Aligner):

    def __init__(self):
        Aligner.__init__(self)
        self.method = 'Perfect_Match'

    def align_library(self, library, template, \
        variant_sequence_quality_threshold = 0, \
        mismatch_quality_threshold = 0,
        output_frequency = 1e5):

        self.template_sequence = template.sequence
        self.variant_sequence_quality_threshold = \
            variant_sequence_quality_threshold
        self.mismatch_quality_threshold = mismatch_quality_threshold
        template_mismatches = 0
        variant_quality_mismatches = 0
        variant_nucleotide_mismatches = 0
        size_mismatches = 0

        template_length = len(self.template_sequence)

        extracted_sequences = []
        uuids = []

        num_sequences = 0

        for fastq_file_name in library.fastq_files:

            fastq_file_name = ws.get_raw_data_path(fastq_file_name)
            fastq_file = open(fastq_file_name, 'r').read().splitlines()

            line_count = 0

            for line in fastq_file:

                if line_count > 10000:
                    break

                if line_count % 4 == 1:
                    if (line_count / 4) % output_frequency == 0:
                        print("Aligned " + str(line_count / 4) + \
                            " sequence(s)")
                    sequence = line

                elif line_count % 4 == 3:

                    quality_string = line
                    
                    if (len(sequence)) == template_length:

                        extracted_sequence, uuid, error_type = \
                            self.extract_sequence(sequence, quality_string)
                        
                        if error_type == 0:
                            extracted_sequences.append(extracted_sequence)
                            uuids.append(uuid)
                        elif error_type == 1:
                            template_mismatches += 1
                        elif error_type == 2:
                            variant_quality_mismatches += 1
                        elif error_type == 3:
                            variant_nucleotide_mismatches += 1
                    else:
                        size_mismatches +=1

                    num_sequences += 1

                line_count += 1

        mismatched_sequences = template_mismatches + variant_quality_mismatches + variant_nucleotide_mismatches + size_mismatches
        
        statistics = {}
        
        statistics["Number of Sequences"] = num_sequences
        statistics["Perfect Match Failure Rate"] = float(mismatched_sequences)/float(num_sequences)
        statistics["Template Mismatch Failure Rate"] = float(template_mismatches)/float(num_sequences)
        statistics["Variant Quality Failure Rate"] = float(variant_quality_mismatches)/float(num_sequences)
        statistics["Variant Nucleotide Mismatch Rate"] = float(variant_nucleotide_mismatches)/float(num_sequences)
        statistics["Template Size Mismatch Rate"] = float(size_mismatches)/float(num_sequences)

        return extracted_sequences, uuids, statistics

    def extract_sequence(self, fastq_line, quality_string):

        extracted_sequence = []
        uuid = []
        num_errors = 0

        template_length = len(self.template_sequence)

        quality_score = FASTQ.convert_quality_string_to_quality_score(quality_string)

        for template_idx in range(template_length):

            # Case 1: The template is an 'I', so this is part of the UUID
            if self.template_sequence[template_idx] == 'I':
                uuid.append(fastq_line[template_idx])

            # Case 2: the template is an N, so this is part of our variant sequence
            elif self.template_sequence[template_idx] == 'N':

                if quality_score[template_idx] < self.variant_sequence_quality_threshold:
                    return False, False, 2

                extracted_sequence.append(fastq_line[template_idx])

            # Case 3: the template is a K, so the variant should be a T or G
            elif self.template_sequence[template_idx] == 'K':

                if quality_score[template_idx] < self.variant_sequence_quality_threshold:
                    return False, False, 2

                if fastq_line[template_idx] != 'G' and fastq_line[template_idx] != 'T':
                    return False, False, 3

                extracted_sequence.append(fastq_line[template_idx])

            # Case 4: the template is an X, so we ignore this
            elif self.template_sequence[template_idx] == 'X':
                # Do nothing
                continue
            # Case 5: This is not part of the variant sequence, so we compare the 
            # template if and only if the quality score is above our threshold
            elif self.template_sequence[template_idx].upper() != fastq_line[template_idx].upper():
                if quality_score[template_idx] > self.mismatch_quality_threshold:
                    return False, False, 1

        return ''.join(extracted_sequence), ''.join(uuid), 0