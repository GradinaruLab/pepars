from workspace import Workspace as ws
from workspace import Database as db
import numpy as np

class Aligner(object):

    def __init__(self):
        pass

    def align(self, alignment, library, progress_callback):

        self._progress_callback = progress_callback

        template_id = alignment.library_templates[library.id]
        template = db.get_template_by_id(template_id)

        if ws.alignment_exists(library, alignment):
            return

        line_count = 0

        for fastq_file_name in library.fastq_files:
            fastq_file_name = ws.get_raw_data_path(fastq_file_name)
            fastq_file = open(fastq_file_name, 'r')

            for line in fastq_file:
                line_count += 1

            fastq_file.close()

        self._num_sequences = line_count / 4

        sequences, uuids, statistics = self.align_library(library, \
            template, **alignment.parameters)

        if len(sequences) != len(uuids) and len(uuids) != 0:
            raise Exception('Number of sequences must match number of \
                UUIDs!')

        sequence_index = 0

        sequence_uuid_counts = {}

        if len(uuids) > 0:
            for sequence_index in range(0, len(sequences)):
                if (sequences[sequence_index], uuids[sequence_index]) not \
                    in sequence_uuid_counts:

                    sequence_uuid_counts[(sequences[sequence_index], \
                        uuids[sequence_index])] = 0

                sequence_uuid_counts[(sequences[sequence_index], \
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

            sequence_uuid_counts_array[sequence_index] = [sequence_uuids[0], sequence_uuids[1], counts]

            sequence_index += 1

        ws.write_sequence_file(library, alignment, \
            sequence_uuid_counts_array)

        alignment.add_statistics(library, statistics)

    def update_num_sequences_aligned(self, num_sequences):

        percent = num_sequences * 100.0 / self._num_sequences

        progress_string = '{0:.2f}'.format(percent) + '% (' + \
            str(num_sequences) + '/' + str(self._num_sequences) + ')'

        self._progress_callback(progress_string)

    @staticmethod
    def validate_alignment(alignment):

        previous_num_variants = -1

        for library, template_id in alignment.library_templates.items():

            template = db.get_template_by_id(template_id)
            num_variants = Aligner.get_num_variants(template)

            if previous_num_variants != -1 and \
                num_variants != previous_num_variants:

                raise Exception('Num variants does not match between all \
                    templates!')

            previous_num_variants = num_variants

    @staticmethod
    def get_num_variants(template):

        num_variants = 0

        variant_markers = set({'N', 'K'})

        for variant_marker in variant_markers:
            num_variants += template.sequence.count(variant_marker)

        return num_variants