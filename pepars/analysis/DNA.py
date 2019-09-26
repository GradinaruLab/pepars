import operator
import numpy
import itertools
import pandas

from ..utils import DNA
from ..utils import Sequence_Trie


def get_template_from_nucleotide_distribution(
        nucleotide_percentiles,
        percentile_threshold=0.05):
    """
    Given an array of nucleotide distribution, return the most likely template
    based on IUPAC grammar.

    :param nucleotide_percentiles: A dictionary of arrays, where each entry
        represents one nucleotide, and the values in the array represent the
        counts or percentiles of that nucleotide in that position
    :param percentile_threshold: What percent a nucleotide must be present to
        be considered part of the template
    :return: A string of length N representing the estimated template
    """

    using_percentiles = False

    for nucleotide, counts in nucleotide_percentiles.items():

        if 0 < counts[0] < 1:
            using_percentiles = True

    template = ""

    if not using_percentiles:

        nucleotide_counts = nucleotide_percentiles

        num_reads = sum([nucleotide_counts[nucleotide][0] for
                         nucleotide in nucleotide_counts])

        nucleotide_percentiles = {}

        for nucleotide in nucleotide_counts:

            nucleotide_percentiles[nucleotide] = []

            for count in nucleotide_counts[nucleotide]:
                nucleotide_percentiles[nucleotide].append(count/num_reads)

    sequence_length = len(nucleotide_percentiles[
                              list(nucleotide_percentiles.keys())[0]])

    for sequence_index in range(sequence_length):

        possible_nucleotides = ""

        for nucleotide in sorted(nucleotide_percentiles):

            if nucleotide_percentiles[nucleotide][sequence_index] > \
                    percentile_threshold:
                possible_nucleotides += nucleotide

        for code, options in DNA.IUPAC_GRAMMAR_MAP.items():
            if options == possible_nucleotides:
                template += code
                break

    return template


def collapse_similar_sequences(sequence_counts,
                               num_nucleotides_off=1):
    """
    Given a list of sequences and their count, recursively collapse ones that
    are similar into their parent(s).

    :param sequence_counts: A list of tuples of sequence and count
    :param num_nucleotides_off: How many nucleotides off a sequence needs to be
        to be considered similar
    :return: A list of tuples of sequence and count
    """

    if num_nucleotides_off != 1:
        raise NotImplementedError("Only support single nucleotide errors")

    print("Sorting sequence counts")
    # Sort list in reverse order (smallest counts first)
    sequence_counts.sort(key=operator.itemgetter(1), reverse=True)

    print("Constructing sequence trie")
    sequence_trie = Sequence_Trie(by_nucleotide=True, allow_invalid=True)

    sequence_index = 0
    for sequence, count in sequence_counts:
        sequence_trie.add(sequence, count)
        sequence_index += 1
        if sequence_index % 100000 == 0:
            print("Added %i sequences" % sequence_index)

    alphabet = set(DNA.get_nucleotides())
    alphabet.add("N")

    for sequence_index, (sequence, count) in enumerate(sequence_counts):

        if sequence_index > 0 and sequence_index % 10000 == 0:
            print("Analyzing sequence %i" % sequence_index)

        min_parent_count = count * 2 - 1

        parent_counts = []

        # Find all possible parent sequences
        for index, character in enumerate(sequence):

            prefix = sequence[0:index]

            parent_node = sequence_trie.get_node(prefix)

            if parent_node is None:
                continue

            for other_character in alphabet.difference(character):

                postfix = other_character + sequence[index + 1:]

                parent_count = parent_node.get_value(postfix)

                if parent_count is None:
                    continue
                elif parent_count < min_parent_count:
                    continue

                parent_counts.append((prefix + postfix, parent_count))

        if len(parent_counts) > 0:

            if len(parent_counts) > 100:
                print("Multiple parents of '%s'" % sequence)
                print(parent_counts)

            parent_count_sum = 0

            for parent, parent_count in parent_counts:
                parent_count_sum += parent_count

            for parent, parent_count in parent_counts:
                count_to_add = parent_count / parent_count_sum * count
                sequence_trie.add(parent, parent_count + count_to_add)
            sequence_trie.add(sequence, 0)

    collapsed_sequence_counts = []

    for sequence, count in sequence_counts:
        new_sequence_count = sequence_trie.get_value(sequence)
        if new_sequence_count > 0:
            collapsed_sequence_counts.append(
                (sequence, int(round(new_sequence_count))))

    return collapsed_sequence_counts


def get_all_sequences_of_distance_n(sequence, n=1):
    """
    Given a sequence, returns a list of all possible sequences that are distance
    n or less away from the given sequence.

    :param sequence: The parent sequence
    :param n: How many nucleotides different each sequence can be
    :return: A list of sequences
    """

    if n != 1:
        raise NotImplementedError("Only doing n=1 for now")

    alphabet = set(DNA.NUCLEOTIDES)

    sequences = []

    # Find all possible parent sequences
    for index, character in enumerate(sequence):

        prefix = sequence[0:index]

        for other_character in alphabet.difference(character):
            postfix = other_character + sequence[index + 1:]

            sequences.append(prefix + postfix)

    return sequences


def check_if_sequence_exists(sequence, sequence_trie, n=1,
                             allow_invalid=True):
    """
    Given a sequence and corresponding sequence trie, check whether any
    sequences exist that are distance n or less away

    :param sequence: The sequence to search for
    :param sequence_trie: The trie to search in
    :param n: How far a sequence can be to be returned
    :param allow_invalid: Whether to consider sequences with Ns
    :return: A list of sequences
    """

    if n != 1 and n != 2:
        raise NotImplementedError("Only doing n={1,2} for now")

    alphabet = set(DNA.get_nucleotides())

    if allow_invalid:
        alphabet.add("N")

    # Find all possible parent sequences
    for index, character in enumerate(sequence):

        prefix = sequence[0:index]

        parent_node = sequence_trie.get_node(prefix)

        if parent_node is None:
            continue

        for other_character in alphabet.difference(character):

            postfix = other_character + sequence[index + 1:]

            exists = parent_node.find(postfix)

            if exists:
                return True

            if n > 1:
                for index_2 in range(index + 1, len(sequence)):

                    midfix = other_character + sequence[index + 1:index_2]

                    middle_node = parent_node.get_child(midfix)

                    if middle_node is None:
                        continue

                    character = sequence[index_2]

                    for other_other_character in alphabet.difference(character):

                        postfix = other_other_character + sequence[index_2 + 1:]

                        exists = middle_node.find(postfix)

                        if exists:
                            return True
    return False


def find_all_sequences_of_distance_n(sequence, sequence_trie, n=1,
                                     allow_invalid=True):
    """
    Given a sequence and corresponding sequence trie, find all sequences
    that are distance n or less away.

    :param sequence: The sequence to search for
    :param sequence_trie: The trie to search in
    :param n: How far a sequence can be to be returned
    :param allow_invalid: Whether to consider sequences with Ns
    :return: A list of sequences
    """

    if n != 1 and n != 2:
        raise NotImplementedError("Only doing n={1,2} for now")

    alphabet = set(DNA.get_nucleotides())

    if allow_invalid:
        alphabet.add("N")

    sequences = []

    # Find all possible parent sequences
    for index, character in enumerate(sequence):

        prefix = sequence[0:index]

        parent_node = sequence_trie.get_node(prefix)

        if parent_node is None:
            continue

        for other_character in alphabet.difference(character):

            postfix = other_character + sequence[index + 1:]

            exists = parent_node.find(postfix)

            if exists:
                sequences.append(prefix + postfix)

            if n > 1:
                for index_2 in range(index + 1, len(sequence)):

                    midfix = other_character + sequence[index + 1:index_2]

                    middle_node = parent_node.get_child(midfix)

                    if middle_node is None:
                        continue

                    character = sequence[index_2]

                    for other_other_character in alphabet.difference(character):

                        postfix = other_other_character + sequence[index_2 + 1:]

                        exists = middle_node.find(postfix)

                        if exists:
                            sequences.append(prefix + midfix + postfix)
    return sequences


def get_amino_acid_probabilities_from_template(template,
                                               allow_stop_codon=False):
    """
    Given a template, return an NxM matrix representing the probability of each
    amino acid in each position

    :param template: A sequence of codons of length N*3

    :param allow_stop_codon: Whether to allow/account for stop codons (True)
        or ignore them (False)

    :return: An NxM numpy array, M is either the number of amino acids (20) or
        the number of amino acids + 1 (for stop codon), if allow_stop_codon
    """

    if len(template) % 3 != 0:
        raise ValueError("Template must be a multiple of 3 in length")

    amino_acids = DNA.get_amino_acids()

    if allow_stop_codon:
        amino_acids.append("#")

    num_amino_acids = len(amino_acids)

    amino_acid_counts = numpy.zeros((int(len(template) / 3), num_amino_acids))

    for i in range(0, len(template), 3):

        possible_nucleotides = [DNA.IUPAC_GRAMMAR_MAP[nucleotide] for
                                nucleotide in template[i:i + 3]]
        for combination in itertools.product(*possible_nucleotides):
            codon = "".join(combination)
            amino_acid = DNA.CODON_AA_MAP[codon]

            if amino_acid == "#":
                if allow_stop_codon:
                    amino_acid_index = -1
                else:
                    continue
            else:
                amino_acid_index = DNA.AMINO_ACID_INDEX_MAP[amino_acid]

            amino_acid_counts[int(i / 3), amino_acid_index] += 1

    amino_acid_probabilities = \
        numpy.divide(amino_acid_counts, amino_acid_counts.sum(axis=1)[:, None])

    amino_acid_probabilities = pandas.DataFrame(
        amino_acid_probabilities.transpose(),
        index=amino_acids,
        columns=[i + 1 for i in range(amino_acid_counts.shape[0])]
    )

    return amino_acid_probabilities


def get_amino_acid_probabilities_from_sequence_counts(
        sequence_counts,
        allow_stop_codon=False):
    """
    Given a dictionary of sequences and their count, returns an NxM matrix
    representing the probability of each amino acid in each position

    :param sequence_counts: A dictionary of sequences and their counts

    :param allow_stop_codon: Whether to allow/account for stop codons (True)
        or ignore them (False)

    :return: An NxM numpy array, M is either the number of amino acids (20) or
        the number of amino acids + 1 (for stop codon), if allow_stop_codon
    """

    amino_acids = DNA.get_amino_acids()

    if allow_stop_codon:
        amino_acids.append("#")

    amino_acid_position_counts = numpy.zeros(
        (len(next(iter(sequence_counts))), len(amino_acids)))

    for sequence, count in sequence_counts.items():

        for character_index, amino_acid in enumerate(sequence):

            if not allow_stop_codon and amino_acid == "#":
                continue

            if amino_acid == "#":
                amino_acid_index = len(amino_acids) - 1
            else:
                amino_acid_index = DNA.AMINO_ACID_INDEX_MAP[amino_acid]

            amino_acid_position_counts[
                character_index, amino_acid_index] += count

    amino_acid_probabilities = \
        amino_acid_position_counts / \
        amino_acid_position_counts.sum(axis=1)[:, None]

    return amino_acid_probabilities


def get_amino_acid_probabilities_from_sequences(
        sequences,
        allow_stop_codon=False):

    amino_acids = DNA.get_amino_acids()

    if allow_stop_codon:
        amino_acids.append("#")

    amino_acid_position_counts = numpy.zeros(
        (len(sequences[0]), len(amino_acids)))

    for sequence in sequences:

        for character_index, amino_acid in enumerate(sequence):

            if not allow_stop_codon and amino_acid == "#":
                continue

            if amino_acid == "#":
                amino_acid_index = len(amino_acids) - 1
            else:
                amino_acid_index = DNA.AMINO_ACID_INDEX_MAP[amino_acid]

            amino_acid_position_counts[
                character_index, amino_acid_index] += 1

    amino_acid_probabilities = \
        amino_acid_position_counts / \
        amino_acid_position_counts.sum(axis=1)[:, None]

    return amino_acid_probabilities
