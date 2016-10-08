from utils import DNA
import itertools

def get_coverage(analysis_set, by_amino_acid = True, use_multiple_templates = False):

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

    sequence_libraries = analysis_set.get_libraries()

    included_sequences = set()
    excluded_sequences = set()

    sequence_library_counts = []

    for sequence_library_name, sequence_library in sequence_libraries.items():
        sequence_library_counts.append(sequence_library.get_sequence_counts(\
            by_amino_acid = by_amino_acid, count_threshold = 0))

    template_index = 0

    iterator = 0

    num_templates = 1

    # Hardcoded for now...
    if use_multiple_templates:
        num_templates = 5

    for template_index in range(0, num_templates):    

        for possible_sequence in itertools.imap(''.join, itertools.product(\
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
                        sequence_to_search_for = possible_sequence + DNA.translate_dna_single('TTGGCGGTGCCTTTTAAGGCACAG')
                    else:
                        sequence_to_search_for = possible_sequence + 'TTGGCGGTGCCTTTTAAGGCACAG'
                elif template_index == 1:
                    if by_amino_acid:
                        sequence_to_search_for = DNA.translate_dna_single('GCCCAA') + possible_sequence + DNA.translate_dna_single('GTGCCTTTTAAGGCACAG')
                    else:
                        sequence_to_search_for = 'GCCCAA' + possible_sequence + 'GTGCCTTTTAAGGCACAG'
                elif template_index == 2:
                    if by_amino_acid:
                        sequence_to_search_for = DNA.translate_dna_single('GCCCAAACTTTG') + possible_sequence + DNA.translate_dna_single('TTTAAGGCACAG')
                    else:
                        sequence_to_search_for = 'GCCCAAACTTTG' + possible_sequence + 'TTTAAGGCACAG'
                elif template_index == 3:
                    if by_amino_acid:
                        sequence_to_search_for = DNA.translate_dna_single('GCCCAAACTTTGGCGGTG') + possible_sequence + DNA.translate_dna_single('GCACAG')
                    else:
                        sequence_to_search_for = 'GCCCAAACTTTGGCGGTG' + possible_sequence + 'GCACAG'
                else:
                    if by_amino_acid:
                        sequence_to_search_for = DNA.translate_dna_single('GCCCAAACTTTGGCGGTGCCTTTT') + possible_sequence
                    else:
                        sequence_to_search_for = 'GCCCAAACTTTGGCGGTGCCTTTT' + possible_sequence

            else:
                sequence_to_search_for = possible_sequence

            sequence_included = False
            for sequence_library_count in sequence_library_counts:
                if sequence_to_search_for in sequence_library_count:
                    included_sequences.add(sequence_to_search_for)
                    sequence_included = True
                    break

            if not sequence_included:
                excluded_sequences.add(sequence_to_search_for)

            iterator += 1

    return included_sequences, excluded_sequences