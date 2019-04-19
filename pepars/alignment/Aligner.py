from ..utils import DNA


class Aligner(object):

    def __init__(self):

        self._num_sequences = 0
        self._progress_callback = None
        self._is_paired_end = False

    def align(self,
              template,
              FASTQ_file_sets,
              alignment_parameters=None,
              progress_callback=None,
              reverse_complement_template=None):

        self._progress_callback = progress_callback

        self._num_sequences = 0

        num_files_per_set = []

        for FASTQ_file_set in FASTQ_file_sets:

            self._num_sequences += FASTQ_file_set.get_read_count()

            num_files_per_set.append(FASTQ_file_set.num_files)

            if FASTQ_file_set.num_files == 1:
                self._is_paired_end = False
            elif FASTQ_file_set.num_files == 2:
                self._is_paired_end = True
            else:
                raise NotImplementedError("Only expecting 1 or 2 files per set")

        if len(set(num_files_per_set)) > 1:
            raise ValueError("All file sets must have the same number of files")

        if alignment_parameters is None:
            alignment_parameters = {}

        sequences, uuids, statistics = self._align(
            FASTQ_file_sets,
            template,
            reverse_complement_template=reverse_complement_template,
            **alignment_parameters)

        if len(sequences) != len(uuids) and len(uuids) != 0:
            raise Exception('Number of sequences must match number of \
                UUIDs!')

        sequence_uuid_counts = {}

        if len(uuids) > 0:
            for sequence_index in range(0, len(sequences)):
                if (sequences[sequence_index], uuids[sequence_index]) not \
                        in sequence_uuid_counts:

                    sequence_uuid_counts[(sequences[sequence_index],
                                          uuids[sequence_index])] = 0

                sequence_uuid_counts[(sequences[sequence_index],
                                      uuids[sequence_index])] += 1
        else:
            for sequence_index in range(0, len(sequences)):

                if (sequences[sequence_index], '') not in \
                        sequence_uuid_counts:

                    sequence_uuid_counts[(sequences[sequence_index], '')] \
                        = 0

                sequence_uuid_counts[(sequences[sequence_index], '')] += 1

        sequence_uuid_counts_array = [[]] * len(sequence_uuid_counts)

        sequence_index = 0

        for sequence_uuids, counts in sequence_uuid_counts.items():

            sequence_uuid_counts_array[sequence_index] = \
                [sequence_uuids[0], sequence_uuids[1], counts]

            sequence_index += 1

        return sequence_uuid_counts_array, statistics

    def update_num_sequences_aligned(self, num_sequences):

        percent = num_sequences * 100.0 / self._num_sequences

        progress_string = "%.2f%% (%i/%i)" % (percent, num_sequences,
                                              self._num_sequences)

        self._progress_callback(progress_string)

    @staticmethod
    def get_num_variants(template):

        num_variants = 0

        for variant_marker in DNA.get_degenerate_nucleotides():
            num_variants += template.sequence.count(variant_marker)

        return num_variants

    def _align(self,
               FASTQ_file_sets,
               template,
               reverse_complement_template,
               alignment_parameters):

        raise NotImplementedError("Aligner subclass must implement _align")
