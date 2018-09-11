import numpy
import math
from utils import DNA
from analysis.Sequence_Library import Sequence_Library
from workspace import Workspace as ws
from . import coverage
import pandas
from scipy import stats
from statsmodels.stats import multitest

def get_probability_of_unseen_sequence(library):

    alignment = ws.get_active_alignment()

    # Check the alignment statistics data to see if we've gotten calculated this before
    if "Unseen Sequence Probability" in alignment.statistics[library.id]:
        return alignment.statistics[library.id]["Unseen Sequence Probability"]

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
    print("num_expected_misreads: %.4f" % num_expected_misreads)
    print("num_single_counts: %i" % num_single_counts)
    print("probability_of_misread_overlap: %.4f" % probability_of_misread_overlap)
    unique_misreads = round(num_expected_misreads * (1-probability_of_misread_overlap))

    n_1 = num_single_counts - unique_misreads

    if n_1 <= num_single_counts * -10:
        raise Exception("Order of magnitude more expected unique misreads than there are single count sequences! This should be impossible")

    if n_1 <= 0:
        n_1 = 1

    probability_unseen = n_1 / alignment.statistics[library.id]["Number of Sequences"]

    alignment.set_statistic(library, "Unseen Sequence Probability", probability_unseen)

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


def find_threshold(labels,
                   top_percent=None,
                   bottom_percent=None,
                   top_count=None,
                   bottom_count=None):

    num_labels = len(labels)

    if top_count is not None:
        top_percent = top_count / num_labels
    elif bottom_count is not None:
        bottom_percent = bottom_count / num_labels

    if top_percent is not None:
        threshold_index = (num_labels - 1) * (1 - top_percent)
        threshold_index = int(math.ceil(threshold_index))
    elif bottom_percent is not None:
        threshold_index = (num_labels - 1) * bottom_percent
        threshold_index = int(math.floor(threshold_index))

    labels_sorted = numpy.sort(labels)

    return labels_sorted[threshold_index]


def get_significance_of_amino_acid_ratios(amino_acid_counts_by_position, amino_acid_biases,
                                          multiple_comparison_correction=True):

    p_values = pandas.DataFrame(numpy.zeros((len(amino_acid_biases.index), len(amino_acid_biases.columns))))
    p_values.index = amino_acid_biases.index
    z_scores = pandas.DataFrame(numpy.zeros((len(amino_acid_biases.index), len(amino_acid_biases.columns))))
    z_scores.index = amino_acid_biases.index

    num_trials = int(amino_acid_counts_by_position.sum(axis=0)[0])

    for amino_acid in amino_acid_biases.index:
        for position_index in range(len(amino_acid_biases.columns)):
            count = int(amino_acid_counts_by_position.loc[amino_acid, position_index])
            p = amino_acid_biases.loc[amino_acid, position_index]
            p_value = stats.binom_test(count, n=num_trials, p=p, alternative="less")
            # z_score = stats.norm.ppf(p_value)
            # z_score = stats.binom.ppf(count, n=num_trials, p=p)
            if p_value == 0:
                z_score = -numpy.inf
            else:
                z_score = -numpy.log10(p_value)
            if p_value > 0.5:
                p_value = stats.binom_test(count, n=num_trials, p=p, alternative="greater")
                # z_score = -stats.norm.ppf(p_value)
                # z_score = -stats.binom.ppf(count, n=num_trials, p=p)
                if p_value == 0:
                    z_score = numpy.inf
                else:
                    z_score = numpy.log10(p_value)
            p_values.loc[amino_acid, position_index] = p_value
            z_scores.loc[amino_acid, position_index] = z_score

    if multiple_comparison_correction:
        _, corrected_p_values, _, _ = multitest.multipletests(
            p_values.values.reshape((p_values.shape[0] * p_values.shape[1],)),
            alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
        corrected_p_values = pandas.DataFrame(corrected_p_values.reshape(p_values.shape))
        corrected_p_values.index = amino_acid_biases.index
        p_values = corrected_p_values

    return p_values, z_scores
