import csv

from ..utils import utils
from ..utils import DNA


def load_sequence_count_file(filename, minimum_num_sequence, group_by_amino_acid):

    sequences = []

    # 1xM list of all the initial counts of all examples
    initial_counts = []

    # MxN list, the counts of each extracted type for each example
    extracted_counts = []

    # Total count of all sequences initially injected into organism
    total_initial_count = 0

    # 1xN list of the sum of all sequences, 1 for each extracted type
    total_extracted_counts = []

    sequences = []
    initial_counts = []
    extracted_counts = []

    with open(filename, 'rb') as sequence_file:
        sequence_reader = csv.reader(sequence_file, delimiter=',')
        headers = next(sequence_reader, None)

        # Loop through all rows and read in their counts
        for row in sequence_reader:

            column_iterator = iter(row)
            sequence = next(column_iterator, None)

            if group_by_amino_acid:
                sequence = DNA.translate_dna_single(sequence)

                # If there are any invalid amino acids, skip them
                if '#' in sequence:
                    continue

            sequences.append(sequence)
            initial_counts.append(int(next(column_iterator, None)))
            extracted_counts.append(int(next(column_iterator, None)))

    sequence_initial_counts = dict()
    sequence_extracted_counts = dict()

    num_sequences = len(sequences)
    print('Num sequences: %d' % num_sequences)

    # Loop through all the sequences we just read in
    for sequence_index in range(0, num_sequences):

        # If this is the first time we're seeing this sequence, initialize our dict
        if sequences[sequence_index] not in sequence_initial_counts:
            sequence_initial_counts[sequences[sequence_index]] = initial_counts[sequence_index]
            sequence_extracted_counts[sequences[sequence_index]] = extracted_counts[sequence_index]
        else:
            sequence_initial_counts[sequences[sequence_index]] += initial_counts[sequence_index]
            sequence_extracted_counts[sequences[sequence_index]] += extracted_counts[sequence_index]

    total_initial_count = 0
    total_extracted_count = 0

    sequences_set = set()

    for sequence in sequence_initial_counts:
        sequences_set.add(sequence)

    # Loop through all our sequences, now that they're nice and grouped uniquely
    for sequence in sequences_set:

        if sequence_initial_counts[sequence] == 0:
            del sequence_initial_counts[sequence]
            del sequence_extracted_counts[sequence]
        elif sequence_initial_counts[sequence] + sequence_extracted_counts[sequence] < minimum_num_sequence:
            del sequence_initial_counts[sequence]
            del sequence_extracted_counts[sequence]
        else:
            total_initial_count += sequence_initial_counts[sequence]
            total_extracted_count += sequence_extracted_counts[sequence]

    sequence_matrix = []

    for sequence in sequence_initial_counts:
        sequence_matrix.append(utils.convert_string_to_char_array(sequence))

    num_sequences = len(sequences)

    # MxN
    fold_enrichments = []

    for sequence in sequence_initial_counts:

        fold_enrichment = (sequence_extracted_counts[sequence]* 1.0 / total_extracted_count) / (sequence_initial_counts[sequence] * 1.0/ total_initial_count)
        fold_enrichments.append([fold_enrichment])

    return sequence_matrix, fold_enrichments


def get_sequences_from_file(file_path):

    sequences = []

    with open(file_path) as sequence_list_file:
        for line in sequence_list_file.readlines():
            sequences.append(line.strip())

    return sequences


def write_sequences_to_file(sequences, file_path):

    with open(file_path, "w") as sequence_list_file:
        for sequence in sequences:
            sequence_list_file.write("%s\n" % sequence)
