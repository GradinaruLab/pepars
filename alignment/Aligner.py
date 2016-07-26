from workspace import Workspace as ws
from workspace import Database as db

class Aligner(object):

    def __init__(self):
        pass

    def align(self, alignment):

        Aligner.validate_alignment(alignment)

        for library_id, template_id in alignment.library_templates.items():

            library = db.get_library_by_id(library_id)
            template = db.get_template_by_id(template_id)

            sequences, uuids, statistics = self.align_library(library, \
                template, **alignment.parameters)

            if len(sequences) != len(uuids) and len(uuids) != 0:
                raise Exception('Number of sequences must match number of \
                    UUIDs!')

            sequence_index = 0

            sequence_uuids = [[0, 0]] * len(sequences)

            if len(uuids) > 0:
                for sequence_index in range(0, len(sequences)):
                    sequence_uuids[sequence_index] = \
                        [sequences[sequence_index], uuids[sequence_index]]
            else:
                for sequence_index in range(0, len(sequences)):
                    sequence_uuids[sequence_index] = \
                        [sequences[sequence_index], []]

            ws.write_sequence_file(library, alignment, sequence_uuids)

            alignment.add_statistics(library, statistics)

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