from utils import DNA
import itertools
from workspace import Workspace as ws
from workspace import Database as db

def get_possible_sequences_by_template(template, by_amino_acid = True):

    if by_amino_acid and len(template) % 3 != 0:
        raise Exception('If using amino acids, \
            template must be a multiple of 3')

    possible_sequence_elements = []
    possible_sequences = set()

    for template_element in template:
        possible_sequence_elements.append(list(DNA.IUPAC[template_element]))

    for possible_sequence in itertools.product(*possible_sequence_elements):
        possible_sequence = ''.join(possible_sequence)

        if by_amino_acid:
            possible_sequence = DNA.translate_dna_single(possible_sequence)

        possible_sequences.add(possible_sequence)

        print(possible_sequence)

    return possible_sequences

def get_possible_sequences(analysis_set, by_amino_acid = True, use_multiple_templates = False):

    sequence_length = analysis_set.get_sequence_length(by_amino_acid)

    # Hardcoded for now
    if use_multiple_templates == True:
        sequence_length = 9

        if by_amino_acid:
            sequence_length /= 3

    if sequence_length == 0:
        return set(), set()

    if by_amino_acid:
        possible_sequence_elements = DNA.get_amino_acids()
    else:
        possible_sequence_elements = DNA.get_nucleotides()


    sequences = set()



    template_index = 0

    num_templates = 1

    # Hardcoded for now...
    if use_multiple_templates:
        num_templates = 5

    for template_index in range(0, num_templates):

        for possible_sequence in map(''.join, itertools.product(\
            possible_sequence_elements,repeat=sequence_length)):

            if not by_amino_acid:

                is_invalid_sequence = False

                for sequence_index in range(0, sequence_length):
                    if sequence_index % 3 == 2 and (possible_sequence[sequence_index] == 'A' or possible_sequence[sequence_index] == 'C'):
                        is_invalid_sequence = True
                        break

                if is_invalid_sequence:
                    continue

            if use_multiple_templates:
                if template_index == 0:
                    if by_amino_acid:
                        sequence = possible_sequence + DNA.translate_dna_single('TTGGCGGTGCCTTTTAAGGCACAG')
                    else:
                        sequence = possible_sequence + 'TTGGCGGTGCCTTTTAAGGCACAG'
                elif template_index == 1:
                    if by_amino_acid:
                        sequence = DNA.translate_dna_single('GCCCAA') + possible_sequence + DNA.translate_dna_single('GTGCCTTTTAAGGCACAG')
                    else:
                        sequence = 'GCCCAA' + possible_sequence + 'GTGCCTTTTAAGGCACAG'
                elif template_index == 2:
                    if by_amino_acid:
                        sequence = DNA.translate_dna_single('GCCCAAACTTTG') + possible_sequence + DNA.translate_dna_single('TTTAAGGCACAG')
                    else:
                        sequence = 'GCCCAAACTTTG' + possible_sequence + 'TTTAAGGCACAG'
                elif template_index == 3:
                    if by_amino_acid:
                        sequence = DNA.translate_dna_single('GCCCAAACTTTGGCGGTG') + possible_sequence + DNA.translate_dna_single('GCACAG')
                    else:
                        sequence = 'GCCCAAACTTTGGCGGTG' + possible_sequence + 'GCACAG'
                else:
                    if by_amino_acid:
                        sequence = DNA.translate_dna_single('GCCCAAACTTTGGCGGTGCCTTTT') + possible_sequence
                    else:
                        sequence = 'GCCCAAACTTTGGCGGTGCCTTTT' + possible_sequence

            else:
                sequence = possible_sequence

            sequences.add(sequence)
    return sequences

def get_coverage_sequences(analysis_set, by_amino_acid = True, use_multiple_templates = False):

    possible_sequences = get_possible_sequences(analysis_set, by_amino_acid, use_multiple_templates)

    included_sequences = set()
    excluded_sequences = set()
    sequence_libraries = analysis_set.get_libraries()

    sequence_library_counts = []

    for sequence_library_name, sequence_library in sequence_libraries.items():
        sequence_library_counts.append(sequence_library.get_sequence_counts(\
            by_amino_acid = by_amino_acid, count_threshold = 0))
    for possible_sequence in possible_sequences:

        sequence_included = False
        for sequence_library_count in sequence_library_counts:
            if possible_sequence in sequence_library_count:
                included_sequences.add(possible_sequence)
                sequence_included = True
                break

        if not sequence_included:
            excluded_sequences.add(possible_sequence)

    return included_sequences, excluded_sequences

def get_coverage(analysis_set, by_amino_acid = True):

    num_included_sequences = 0

    variant_template = ''

    alignment = ws.get_active_alignment()

    # Get the templates for each analysis set
    for library_name in analysis_set.get_libraries():
        db_library = db.get_library(library_name)
        template_id = alignment.library_templates[db_library.id]
        template = db.get_template_by_id(template_id)
        if len(variant_template) == 0:
            variant_template = template.get_variant_template()
        elif variant_template != template.get_variant_template():
            raise Exception('Variant templates must match in an analysis set to do coverage analysis!')

    num_possible_sequences = 1

    if by_amino_acid:
        if len(variant_template) % 3 != 0:
            raise Exception('Can\'t analyze by amino acid when variant sequence isn\'t groups of 3!')
        variant_length = int(len(variant_template) / 3)

        # For now, assume all amino acids are possible. Should do something with IUPAC later
        for variant_index in range(0, variant_length):
            num_possible_sequences *= 20;
    else:
        for variant_index in range(0, len(variant_template)):
            num_possible_sequences *= len(DNA.IUPAC[variant_template[variant_index]])

    # We assume all sequences in the analysis set have been aligned against the template,
    # so they must match the template. So, all unique sequences are the included sequences

    included_sequences = set()

    for library_name, library in analysis_set.get_libraries().items():
        sequence_counts = library.get_sequence_counts(by_amino_acid, count_threshold = 0)

        for sequence, count in sequence_counts.items():
            included_sequences.add(sequence)

    return len(included_sequences), num_possible_sequences
