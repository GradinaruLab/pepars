from enum import Enum


class Confidence_Metric(Enum):

    LOG_COUNTS = 1
    LOG_PROBABILITY = 2


def get_sequence_confidence(
        sequence_counts,
        confidence_metric=Confidence_Metric.LOG_COUNTS):
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

