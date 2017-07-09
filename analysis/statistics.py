import numpy
import math
from utils import DNA
from analysis.Sequence_Library import Sequence_Library
from workspace import Workspace as ws
from . import coverage

def get_probability_of_unseen_sequence(library):

    alignment = ws.get_active_alignment()

    # Check the alignment statistics data to see if we've gotten calculated this before
    if "Unseen Sequence Probability" in alignment.statistics[library.id]:
        return alignment.statistics[library.id]["Unseeen Sequence Probability"]

    if "Expected Number of Misreads" not in alignment.statistics[library.id]:
        raise Exception("Missing 'Expected Number of Misreads' from statistics. Align with an alignment method that generates ")

    sequence_library = Sequence_Library(library)

    num_expected_misreads = alignment.statistics[library.id]["Expected Number of Misreads"]
    sequence_counts = sequence_library.get_sequence_counts(by_amino_acid=False, count_threshold=0, filter_invalid=True)
    num_single_counts = 0

    for sequence, count in sequence_counts.items():
        if count == 1:
            num_single_counts += 1

    probability_of_misread_overlap = coverage.get_probability_of_single_misread_existing(library)
    unique_misreads = round(num_expected_misreads * (1-probability_of_misread_overlap))
    num_expected_misreads = round(num_expected_misreads)

    n_1 = num_single_counts - unique_misreads

    if n_1 <= 0:
        raise Exception("More expected unique misreads than there are single count sequences! This should be impossible")

    probability_unseen = n_1 / alignment.statistics[library.id]["Number of Sequences"]

    alignment.set_statistic(library, "Unseeen Sequence Probability", probability_unseen)

    return probability_unseen

def get_sequence_count_bins(sequence_library, bins = None):

    if not bins:
        bins = [1, 5, 9, 49, 99, 499, 999, 4999, 9999]

    sequence_counts = sequence_library.get_sequence_counts(count_threshold = 0)

    return get_bin_counts(list(sequence_counts.values()), bins)

def get_enrichment_bins(sequence_enrichments, bins = None):

    if not bins:
        bins = [-4, -3, -2, -1, 0, 1, 2, 3, 4]

    return get_bin_counts(list(sequence_enrichments.values()), bins)

def get_bin_counts(values, bins):

    bins.append(0)

    bin_counts = [0] * len(bins)

    for value in values:

        bin_index = 0

        while value > bins[bin_index] and bin_index < len(bins) - 1:
            bin_index += 1

        bin_counts[bin_index] += 1

    return bins, bin_counts

def binomial(n, p, k):

    answer = 1
    for i in range(k):
        answer *= (n - i)

    answer /= math.factorial(k)

    answer *= pow(p, k)
    answer *= pow(1-p, n-k)

    return answer

def get_expected_count_distribution(template, num_samples):

    num_possible_sequences = 1

    for character in template:
        num_possible_sequences *= len(list(DNA.IUPAC[character]))

    k = 0

    probability_of_sequence = 1/num_possible_sequences

    probability_unseen = binomial(num_samples, probability_of_sequence, k)

    probability_seen = 1 - probability_unseen

    distribution = {}

    while True:

        k += 1

        probability = binomial(num_samples, probability_of_sequence, k)

        relative_probability = probability / probability_seen

        num_expected_sequences = round(relative_probability * num_samples)

        if num_expected_sequences >= 1:
            distribution[k] = num_expected_sequences
        else:
            break

    return distribution