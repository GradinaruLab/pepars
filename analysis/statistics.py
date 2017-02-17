import numpy
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