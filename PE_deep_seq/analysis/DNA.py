import operator

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

        for nucleotide in nucleotide_percentiles:

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

    # Sort list in reverse order (smallest counts first)
    sequence_counts.sort(key=operator.itemgetter(1), reverse=True)

    print("Building sequence trie")
    sequence_trie = Sequence_Trie(by_nucleotide=True, allow_invalid=True)
    for sequence, count in sequence_counts:
        sequence_trie.add(sequence, count)

    alphabet = set(DNA.get_nucleotides())
    alphabet.add("N")

    for sequence_index, (sequence, count) in enumerate(sequence_counts):

        if sequence_index % 10000 == 0:
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
