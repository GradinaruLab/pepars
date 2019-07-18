from .Aligner import Aligner
from ..utils import FASTQ
from ..utils import DNA


class Perfect_Match_Aligner(Aligner):

    def __init__(self):
        Aligner.__init__(self)
        self.method = self.__class__.__name__
        self._template = None
        self._reverse_complement_template = None
        self._variant_sequence_quality_threshold = None
        self._mismatch_quality_threshold = None

    def _align(self,
               FASTQ_file_sets,
               template,
               reverse_complement_template,
               variant_sequence_quality_threshold=0,
               mismatch_quality_threshold=0,
               output_frequency=1e5):

        self._template = template
        self._reverse_complement_template = reverse_complement_template
        self._variant_sequence_quality_threshold = \
            int(variant_sequence_quality_threshold)
        self._mismatch_quality_threshold = int(mismatch_quality_threshold)
        template_mismatches = 0
        variant_quality_mismatches = 0
        variant_nucleotide_mismatches = 0
        size_mismatches = 0
        invalid_nucleotides = 0
        paired_end_mismatches = 0

        template_length = len(self._template)

        extracted_sequences = []
        uuids = []
        error_probabilities = []
        sequence_counts = {}

        num_sequences = 0
 
        for FASTQ_file_set in FASTQ_file_sets:

            if self._is_paired_end:
                self._progress_callback("Reading '%s'" %
                                        FASTQ_file_set.files[0].file_name)
                self._progress_callback("Reading '%s'" %
                                        FASTQ_file_set.files[1].file_name)

            sequence_count = 0

            FASTQ_sequence_index = 0

            for sequences, quality_strings in \
                    FASTQ_file_set.get_sequence_quality_iterator():

                if sequence_count % output_frequency == 0:
                    self.update_num_sequences_aligned(num_sequences)

                if (len(sequences[0])) == template_length:

                    extracted_sequence, uuid, error_type, error_probability = \
                        self.extract_sequence(sequences[0], quality_strings[0],
                                              self._template)
                else:
                    error_type = 6

                # Only bother checking paired end if error_type is 0
                if error_type == 0 and self._is_paired_end:

                    if len(sequences[1]) != \
                           len(self._reverse_complement_template):
                        error_type = 6
                    else:
                        extracted_sequence_r2, _, error_type,\
                            error_probability_r2 = \
                            self.extract_sequence(
                                sequences[1],
                                quality_strings[1],
                                self._reverse_complement_template)

                        if error_type == 0:
                            extracted_sequence_r2 = \
                                DNA.translate_reverse_complement(
                                    extracted_sequence_r2)
                            if extracted_sequence_r2 != extracted_sequence:
                                error_type = 5

                            error_probability *= error_probability_r2

                if error_type == 0:

                    extracted_sequences.append(extracted_sequence)
                    if extracted_sequence in sequence_counts:
                        sequence_counts[extracted_sequence] += 1
                    else:
                        sequence_counts[extracted_sequence] = 1
                    uuids.append(uuid)
                    error_probabilities.append(error_probability)

                if error_type == 1:
                    template_mismatches += 1
                elif error_type == 2:
                    variant_quality_mismatches += 1
                elif error_type == 3:
                    variant_nucleotide_mismatches += 1
                elif error_type == 4:
                    invalid_nucleotides += 1
                elif error_type == 5:
                    paired_end_mismatches += 1
                elif error_type == 6:
                    size_mismatches += 1

                sequence_count += 1

                FASTQ_sequence_index += 1

                num_sequences += 1

        mismatched_sequences = template_mismatches + variant_quality_mismatches\
            + variant_nucleotide_mismatches + size_mismatches\
            + paired_end_mismatches + invalid_nucleotides

        num_single_counts = 0
        for sequence, count in sequence_counts.items():
            if count == 1:
                num_single_counts += 1
        
        statistics = {}
        
        statistics["Number of Sequences"] = num_sequences
        statistics["Alignment rate"] = float(
            1.0-float(mismatched_sequences)/float(num_sequences))
        statistics["Template Mismatch Failure Rate"] = float(
            template_mismatches)/float(num_sequences)
        statistics["Paired End Mismatch Failure Rate"] = float(
            paired_end_mismatches)/float(num_sequences)
        statistics["Variant Quality Failure Rate"] = float(
            variant_quality_mismatches)/float(num_sequences)
        statistics["Variant Nucleotide Mismatch Rate"] = float(
            variant_nucleotide_mismatches)/float(num_sequences)
        statistics["Template Size Mismatch Rate"] = float(
            size_mismatches)/float(num_sequences)
        statistics["Invalid nucleotide Rate"] = float(
            invalid_nucleotides)/float(num_sequences)
        statistics["Expected Number of Misreads"] = sum(error_probabilities)
        statistics["Number Single Count Sequences"] = num_single_counts
        if self._is_paired_end:
            statistics["Paired end match rate"] = \
                1 - (paired_end_mismatches / num_sequences)

        return extracted_sequences, uuids, statistics

    def extract_sequence(self, fastq_line, quality_string, template_sequence):

        extracted_sequence = []
        uuid = []
        probability_of_perfect_read = 1.0

        template_length = len(template_sequence)

        quality_score = \
            FASTQ.convert_quality_string_to_quality_score(quality_string)

        for template_idx in range(template_length):

            # Case 1: the template is an X, so we ignore this
            if template_sequence[template_idx] == 'X':
                # Do nothing
                continue
            # Case 2: The template is looking for a specific nudleotide, so we
            # compare the template if and only if the quality score is above our
            # threshold
            elif template_sequence[template_idx] in DNA.get_nucleotides():

                # If the quality score is below our threshold, we ignore this
                # comparison
                if quality_score[template_idx] < \
                        self._mismatch_quality_threshold:
                    continue

                # If this isn't a match, it's an error
                if template_sequence[template_idx].upper() != \
                        fastq_line[template_idx].upper():
                    return False, False, 1, 1
            # Case 3: If this is not a valid nucleotide, it's an error
            elif fastq_line[template_idx] not in DNA.get_nucleotides():
                return False, False, 4, 1
            # Case 4: The template is an 'I', so this is part of the UUID
            elif template_sequence[template_idx] == 'I':
                uuid.append(fastq_line[template_idx])

            # Case 5: the template is an N, so this is part of our variant
            # sequence
            else:
                if template_sequence[template_idx] == 'N':

                    if quality_score[template_idx] < \
                            int(self._variant_sequence_quality_threshold):
                        return False, False, 2, 1

                # Case 6: the template is a K, so the variant should be a T or G
                elif template_sequence[template_idx] == 'K':

                    if quality_score[template_idx] < \
                            self._variant_sequence_quality_threshold:
                        return False, False, 2, 1

                    if fastq_line[template_idx] != 'G' and \
                            fastq_line[template_idx] != 'T':
                        return False, False, 3, 1

                probability_of_error = \
                    FASTQ.convert_quality_score_to_probability_of_error(
                        quality_score[template_idx])

                extracted_sequence.append(fastq_line[template_idx])
                probability_of_perfect_read *= (1.0 - probability_of_error)

        return ''.join(extracted_sequence),\
               ''.join(uuid),\
               0,\
               1.0 - probability_of_perfect_read
