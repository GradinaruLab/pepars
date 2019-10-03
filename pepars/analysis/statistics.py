import numpy
import math
from ..utils import DNA
import pandas
from statsmodels.stats import multitest
from statsmodels.stats import proportion

from enum import Enum


class Test_Type(Enum):

    BINOMIAL_NORMAL_APPROXIMATION_LOG2 = 1
    BINOMIAL_NORMAL_APPROXIMATION_Z = 2


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


def get_significance_of_amino_acid_ratios(
        amino_acid_counts_by_position_output,
        amino_acid_counts_by_position_input,
        multiple_comparison_correction=True,
        test_type=Test_Type.BINOMIAL_NORMAL_APPROXIMATION_LOG2
):
    """

    :param amino_acid_counts_by_position_output: a pandas DataFrame of amino
        acid counts, where rows are the amino acids, and columns are the
        positions in the sequence
    :param amino_acid_counts_by_position_input: The expected amino acid
        probabilities. A pandas DataFrame, where rows are the amino acids, and
        columns are the positions in the sequence
    :param multiple_comparison_correction: Whether to perform multiple
        comparison correction
    :param test_type: The type of statistical test to perform. See
        statistics.Test_Type
    :return: (p_values, magnitude_scores)
        The p_values at each amino acid/position and their associated
        magnitude scores. The meaning of this score is dependent on test_type
    """

    p_values = pandas.DataFrame(
        numpy.zeros((len(amino_acid_counts_by_position_input.index),
                     len(amino_acid_counts_by_position_input.columns))),
        index=amino_acid_counts_by_position_input.index,
        columns=amino_acid_counts_by_position_input.columns
    )

    magnitude_scores = pandas.DataFrame(
        numpy.zeros((len(amino_acid_counts_by_position_input.index),
                     len(amino_acid_counts_by_position_input.columns))),
        index=amino_acid_counts_by_position_input.index,
        columns=amino_acid_counts_by_position_input.columns
    )

    num_trials_output = amino_acid_counts_by_position_output.sum(axis=0)[1]
    num_trials_input = amino_acid_counts_by_position_input.sum(axis=0)[1]

    amino_acid_biases = amino_acid_counts_by_position_input /\
        amino_acid_counts_by_position_input.sum(axis=0)

    # Check if the user has provided probabilities instead of counts
    if round(num_trials_input) == 1.0:
        is_single_sample = True
    else:
        is_single_sample = False

    for amino_acid in p_values.index:
        for position_index in p_values.columns:

            if is_single_sample:
                amino_acid_count = round(amino_acid_counts_by_position_output.loc[
                                           amino_acid, position_index])

                if amino_acid_count == 0:
                    z_score = numpy.nan
                    p_value = numpy.nan
                else:
                    amino_acid_bias = amino_acid_biases.loc[
                        amino_acid, position_index]
                    z_score, p_value = proportion.proportions_ztest(
                        amino_acid_count,
                        num_trials_output,
                        value=amino_acid_bias
                    )
            else:
                amino_acid_count_output = int(
                    amino_acid_counts_by_position_output.loc[
                        amino_acid, position_index])
                amino_acid_count_input = int(
                    amino_acid_counts_by_position_input.loc[
                        amino_acid, position_index])

                if amino_acid_count_input == 0 and \
                        amino_acid_count_output == 0:
                    z_score = numpy.nan
                    p_value = numpy.nan
                else:
                    z_score, p_value = proportion.proportions_ztest(
                        [amino_acid_count_output, amino_acid_count_input],
                        [num_trials_output, num_trials_input],
                        value=0
                    )

            p_values.loc[amino_acid, position_index] = p_value
            magnitude_scores.loc[amino_acid, position_index] = z_score

    if test_type == Test_Type.BINOMIAL_NORMAL_APPROXIMATION_LOG2:

        amino_acid_ratios = amino_acid_counts_by_position_output /\
                            num_trials_output
        amino_acid_biases[amino_acid_biases == 0] = \
            amino_acid_biases[amino_acid_biases != 0].min().min() / 2
        amino_acid_ratios[amino_acid_ratios == 0] = \
            amino_acid_ratios[amino_acid_ratios != 0].min().min() / 2
        magnitude_scores = amino_acid_ratios / amino_acid_biases
        magnitude_scores = numpy.log2(magnitude_scores)

    if multiple_comparison_correction:

        p_values_array =  \
            p_values.values.reshape((p_values.shape[0] * p_values.shape[1],))

        p_values_array_not_nan = p_values_array[~numpy.isnan(p_values_array)]
        _, corrected_p_values, _, _ = multitest.multipletests(
            p_values_array_not_nan, alpha=0.05, method='b', is_sorted=False,
            returnsorted=False)

        p_values_array[~numpy.isnan(p_values_array)] = corrected_p_values

        corrected_p_values = pandas.DataFrame(
            p_values_array.reshape(p_values.shape),
            index=amino_acid_biases.index,
            columns=amino_acid_biases.columns)
        p_values = corrected_p_values

    return p_values, magnitude_scores
