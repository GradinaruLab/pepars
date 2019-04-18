import numpy
import math
from ..utils import DNA
import pandas
from scipy import stats
from statsmodels.stats import multitest
from enum import Enum


class Test_Type(Enum):

    BINOMIAL_NORMAL_APPROXIMATION = 1
    BINOMIAL_LOG_SCORE = 2
    NORMAL_ASSUMPTION = 3
    LOG_FOLD_CHANGE = 4
    STANDARDIZATION = 5

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
                                          multiple_comparison_correction=True,
                                          test_type=Test_Type.BINOMIAL_NORMAL_APPROXIMATION):

    p_values = pandas.DataFrame(
        numpy.zeros((len(amino_acid_biases.index), len(amino_acid_biases.columns))),
        index=amino_acid_biases.index,
        columns=amino_acid_biases.columns)
    z_scores = pandas.DataFrame(numpy.zeros((len(amino_acid_biases.index), len(amino_acid_biases.columns))),
        index=amino_acid_biases.index,
        columns=amino_acid_biases.columns)
    z_scores.index = amino_acid_biases.index
    p_values.columns = amino_acid_biases.columns

    num_trials = int(amino_acid_counts_by_position.sum(axis=0)[1])

    amino_acid_ratios = amino_acid_counts_by_position / num_trials
    amino_acid_ratios = amino_acid_ratios/amino_acid_biases
    amino_acid_ratios_log = numpy.log2(amino_acid_ratios)
    amino_acid_ratios_log[amino_acid_ratios_log == -numpy.inf] = amino_acid_ratios_log[amino_acid_ratios_log != -numpy.inf].min().min() - 1

    for amino_acid in amino_acid_biases.index:
        for position_index in amino_acid_biases.columns:
            count = int(amino_acid_counts_by_position.loc[amino_acid, position_index])
            p = amino_acid_biases.loc[amino_acid, position_index]

            if test_type == Test_Type.BINOMIAL_NORMAL_APPROXIMATION:
                q = 1 - p
                standard_deviation = numpy.sqrt(num_trials * p * q)
                expected_value = num_trials * p
                z_score = (count - expected_value) / standard_deviation
                p_value = stats.norm.cdf(z_score)

                # if p_value > 0.5:
                #     p_value = (1 - p_value)

            elif test_type == Test_Type.BINOMIAL_LOG_SCORE:
                p_value = stats.binom_test(count, n=num_trials, p=p, alternative="less")
                if p_value == 0:
                    z_score = numpy.inf
                else:
                    z_score = numpy.log10(p_value)
                if p_value > 0.5:
                    p_value = stats.binom_test(count, n=num_trials, p=p, alternative="greater")
                    if p_value == 0:
                        z_score = -numpy.inf
                    else:
                        z_score = -numpy.log10(p_value)
            elif test_type == Test_Type.NORMAL_ASSUMPTION:

                ratio = amino_acid_ratios_log.loc[amino_acid, position_index]
                z_score = ratio/amino_acid_ratios_log.values.std()
                p_value = stats.norm.cdf(z_score)

            elif test_type == Test_Type.LOG_FOLD_CHANGE:

                ratio = amino_acid_ratios_log.loc[amino_acid, position_index]
                z_score = ratio/amino_acid_ratios_log.values.std()
                p_value = stats.norm.cdf(z_score)
                z_score = ratio

            elif test_type == Test_Type.STANDARDIZATION:

                z_score = (amino_acid_ratios.loc[amino_acid, position_index] - amino_acid_ratios.values.mean()) / amino_acid_ratios.values.std()
                p_value = stats.norm.cdf(z_score)

            p_values.loc[amino_acid, position_index] = p_value
            z_scores.loc[amino_acid, position_index] = z_score

    z_scores[z_scores == -numpy.inf] = z_scores[z_scores != -numpy.inf].min().min() - 1
    z_scores[z_scores == numpy.inf] = z_scores[z_scores != numpy.inf].max().max() + 1

    if multiple_comparison_correction:
        _, corrected_p_values, _, _ = multitest.multipletests(
            p_values.values.reshape((p_values.shape[0] * p_values.shape[1],)),
            alpha=0.05, method='b', is_sorted=False, returnsorted=False)
        corrected_p_values = pandas.DataFrame(corrected_p_values.reshape(p_values.shape),
            index=amino_acid_biases.index,
            columns=amino_acid_biases.columns)
        p_values = corrected_p_values

    return p_values, z_scores
