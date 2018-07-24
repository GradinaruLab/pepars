from enum import Enum
import numpy

# Log offset is 2 so that we don't get really low/zero values for the
# zero count sequences
LOG_OFFSET = 2


class Confidence_Metric(Enum):

    LOG_COUNTS = 1
    LOG_PROBABILITY = 2
    LINEAR = 3
    LOG_LOG_COUNTS = 4


def get_sequence_confidences(
        sequence_counts,
        confidence_metric=Confidence_Metric.LOG_COUNTS,
        normalize=True):
    """ Returns the confidence of each sequence in a

    Required Arguments
    ------------------
    sequence_counts : numpy.ndarray
        An NxD array of sequence counts, where each row is a sequence, and each
        column is its count in different samples

    Optional Arguments
    -----------------
    confidence_metric : Confidence_Metric
        Which confidence metric to use, see Confidence_Metric for details.
        Default is Confidence_Metric.LOG_COUNTS

    Returns
    -------
    numpy.ndarray
        An Nx1 array of sequence confidences, retains the order of sequences
    """

    if sequence_counts[sequence_counts == 0].sum() > 0:
        raise ValueError("Some sequence counts have a count of 0, "
                         "cannot compute confidence")

    if confidence_metric == Confidence_Metric.LOG_LOG_COUNTS:

        confidences = numpy.zeros((sequence_counts.shape[0],))

        for column_index in range(sequence_counts.shape[1]):

            confidences += numpy.log(
                numpy.log(sequence_counts[:, column_index] + LOG_OFFSET) /
                numpy.log(max(sequence_counts[:, column_index]) + LOG_OFFSET))

    elif confidence_metric == Confidence_Metric.LOG_COUNTS:

        confidences = numpy.ones((sequence_counts.shape[0],))

        for column_index in range(sequence_counts.shape[1]):

            confidences *= \
                numpy.log(sequence_counts[:, column_index] + LOG_OFFSET) / \
                numpy.log(max(sequence_counts[:, column_index]) + LOG_OFFSET)

    elif confidence_metric == Confidence_Metric.LOG_PROBABILITY:

        confidences = numpy.zeros((sequence_counts.shape[0],))

        for column_index in range(sequence_counts.shape[1]):

            confidences += numpy.log(
                (sequence_counts[:, column_index]) /
                max(sequence_counts[:, column_index]))

    elif confidence_metric == Confidence_Metric.LINEAR:

        confidences = numpy.zeros((sequence_counts.shape[0],))

        for column_index in range(sequence_counts.shape[1]):

            confidences += (sequence_counts[:, column_index]) / \
                max(sequence_counts[:, column_index])

    # Normalize it
    if normalize:
        confidences = confidences - confidences.min()
        confidences = confidences / confidences.max()

    return confidences
