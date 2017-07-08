import numpy
import math
from utils import DNA
from analysis.Sequence_Library import Sequence_Library

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

        #print("Probability of seeing a sequence %i time(s): %0.4f%%" % (k, probability * 100))

        relative_probability = probability / probability_seen

        #print("Relative probability of seeing a sequence %i time(s): %0.4f%%" % (k, relative_probability * 100))

        num_expected_sequences = round(relative_probability * num_samples)

        #print("Expected number of sequences that are present %i time(s): %i" % (k, num_expected_sequences))

        if num_expected_sequences >= 1:
            distribution[k] = num_expected_sequences
        else:
            break

    return distribution