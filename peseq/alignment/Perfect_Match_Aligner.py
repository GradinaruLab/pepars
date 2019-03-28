from .Aligner import Aligner
from ..utils import FASTQ
from protfarm.workspace import Workspace as ws
from protfarm.workspace import Database as db
from ..utils import DNA

class Perfect_Match_Aligner(Aligner):

    def __init__(self):
        Aligner.__init__(self)
        self.method = self.__class__.__name__

    def align_library(self, library, template, \
        variant_sequence_quality_threshold = 0, \
        mismatch_quality_threshold = 0,
        output_frequency = 1e5):

        self.template = template
        self.template_sequence = template.sequence
        self.variant_sequence_quality_threshold = \
            int(variant_sequence_quality_threshold)
        self.mismatch_quality_threshold = int(mismatch_quality_threshold)
        template_mismatches = 0
        variant_quality_mismatches = 0
        variant_nucleotide_mismatches = 0
        size_mismatches = 0
        invalid_nucleotides = 0
        paired_end_mismatches = 0

        print('Aligning library \'' + library.name + '\'')

        template_length = len(self.template_sequence)

        extracted_sequences = []
        uuids = []
        error_probabilities = []
        sequence_counts = {}

        num_sequences = 0

        # Find if any files are the reverse complement

        is_paired_end = False

        for fastq_file_name in library.fastq_files:

            FASTQ_file_object = db.get_FASTQ_file(fastq_file_name)

            if FASTQ_file_object.is_reverse_complement:
                is_paired_end = True

        if is_paired_end:
            self._num_sequences = int(self._num_sequences / 2)

        previous_FASTQ_sequences = []
        paired_end_comparisons = 0
 
        for fastq_file_index in range(len(library.fastq_files)):

            fastq_file_name = library.fastq_files[fastq_file_index]

            FASTQ_file_object = db.get_FASTQ_file(fastq_file_name)

            self._progress_callback("Reading '%s'" % fastq_file_name)

            if FASTQ_file_object.is_reverse_complement:
                template = db.get_template_by_id(self.template.reverse_complement_template_id)
                self.template_sequence = template.sequence
            else:
                self.template_sequence = self.template.sequence

            fastq_file = ws.get_fastq_file(fastq_file_name)

            sequence_count = 0

            FASTQ_sequence_index = 0

            for sequence, quality_string in \
                    fastq_file.get_sequence_quality_iterator():

                if sequence_count % output_frequency == 0:
                    self.update_num_sequences_aligned(num_sequences)

                if (len(sequence)) == template_length:

                    extracted_sequence, uuid, error_type, error_probability = \
                        self.extract_sequence(sequence, quality_string)

                    if is_paired_end:

                        if error_type == 0 and FASTQ_file_object.is_reverse_complement:
                            extracted_sequence = DNA.get_reverse_complement(extracted_sequence)

                        # If this is paired end, we just log it and wait for the next pass
                        if fastq_file_index == 0:
                            sequence_error = [extracted_sequence, error_type, error_probability, uuid]
                            previous_FASTQ_sequences.append(sequence_error)
                            # Not really an error, but we want to skip it
                            error_type = 6
                        else:
                            # If the previous FASTQ file had an error, we just use the new one
                            if previous_FASTQ_sequences[FASTQ_sequence_index][1] != 0:
                                pass
                            # If the current FASTQ file has an error, we just use the previous one
                            elif error_type != 0:
                                extracted_sequence = previous_FASTQ_sequences[FASTQ_sequence_index][0]
                                uuid = previous_FASTQ_sequences[FASTQ_sequence_index][3]
                                error_type = previous_FASTQ_sequences[FASTQ_sequence_index][1]
                                error_probability = previous_FASTQ_sequences[FASTQ_sequence_index][2]
                            # If neither had an error, we compare and only include if they match
                            else:
                                paired_end_comparisons += 1
                                if previous_FASTQ_sequences[FASTQ_sequence_index][0] != extracted_sequence:
                                    error_type = 5
                                # If the match, we have even greater confidence in them
                                else:
                                    error_probability = error_probability * previous_FASTQ_sequences[FASTQ_sequence_index][2]

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
                    elif error_type == 5:
                        paired_end_mismatches += 1
                else:
                    size_mismatches +=1

                sequence_count += 1

                FASTQ_sequence_index += 1

                if not is_paired_end or fastq_file_index == 1:
                    num_sequences += 1

            ws.close_fastq_file(fastq_file_name)

        mismatched_sequences = template_mismatches + variant_quality_mismatches\
            + variant_nucleotide_mismatches + size_mismatches + paired_end_mismatches + invalid_nucleotides

        num_single_counts = 0
        for sequence, count in sequence_counts.items():
            if count == 1:
                num_single_counts += 1
        
        statistics = {}
        
        statistics["Number of Sequences"] = num_sequences
        statistics["Alignment rate"] = float(1.0-float(mismatched_sequences)/float(num_sequences))
        statistics["Template Mismatch Failure Rate"] = float(template_mismatches)/float(num_sequences)
        statistics["Paired End Mismatch Failure Rate"] = float(paired_end_mismatches)/float(num_sequences)
        statistics["Variant Quality Failure Rate"] = float(variant_quality_mismatches)/float(num_sequences)
        statistics["Variant Nucleotide Mismatch Rate"] = float(variant_nucleotide_mismatches)/float(num_sequences)
        statistics["Template Size Mismatch Rate"] = float(size_mismatches)/float(num_sequences)
        statistics["Invalid nucleotide Rate"] = float(invalid_nucleotides)/float(num_sequences)
        statistics["Expected Number of Misreads"] = sum(error_probabilities)
        statistics["Number Single Count Sequences"] = num_single_counts
        if paired_end_comparisons > 0:
            statistics["Paired end match rate"] = 1 - (paired_end_mismatches / paired_end_comparisons) 

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