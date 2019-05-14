import numpy

from ..utils import utils
from ..utils import FASTQ_File
from ..utils import FASTQ as FASTQ_utils
from ..utils import DNA as DNA_utils


def get_read_length(FASTQ_file_path):
    """
    Get the read length from a FASTQ file. Assumes all reads are of the same
    length

    :param FASTQ_file_path: Path to a FASTQ file
    :return: The length of reads, as determined by the first read in the file
    """

    file = FASTQ_File(FASTQ_file_path)

    read_length = None

    for sequence in file.get_sequence_iterator():
        read_length = len(sequence)
        break

    file.close()

    return read_length


def get_quality_score_distribution(FASTQ_file_path):
    """
    Given a FASTQ file path, return its distribution of quality scores in each
    position

    :param FASTQ_file_path: Path to a FASTQ file
    :return: An array of dictionaries, one entry for each position in the
        sequencing read. Each dictionary entry is the quality score and the
        values are the counts of that score at that position.
    """

    read_length = get_read_length(FASTQ_file_path)

    file = FASTQ_File(FASTQ_file_path)

    quality_score_counts_dict = [
        {score: 0 for score in range(FASTQ_utils.MAX_QUALITY_SCORE + 1)} for _
        in range(read_length)]

    quality_score_offset = FASTQ_utils.QUALITY_SCORE_OFFSET

    for quality_string in file.get_quality_iterator():

        quality_vector = [ord(quality_score) - quality_score_offset for
                          quality_score in quality_string]

        for position, quality_score in enumerate(quality_vector):
            quality_score_counts_dict[position][quality_score] += 1

    quality_score_counts = numpy.zeros(
        (read_length, FASTQ_utils.MAX_QUALITY_SCORE + 1))

    for position in range(read_length):
        for quality_score, count in quality_score_counts_dict[position].items():
            quality_score_counts[position, quality_score] = count

    return quality_score_counts


def get_nucleotide_distribution(FASTQ_file_path, include_N=True):
    """
    Given a FASTQ file path, return its distribution of nucleotides in each
    position

    :param FASTQ_file_path: Path to a FASTQ file
    :param include_N: Whether to include Ns in the counts
    :return: A dictionary of arrays, one entry for each unique character in the
        sequencing reads. Each array is the length of the reads, and the entries
        are the counts of that character at that position.
    """

    read_length = get_read_length(FASTQ_file_path)

    nucleotides = DNA_utils.get_nucleotides().copy()
    nucleotides.append("N")

    # Initialize the array - nucleotides + 1 for N
    sequence_counts = {character: [0 for _ in range(read_length)] for
                       character in nucleotides}

    file = FASTQ_File(FASTQ_file_path)

    for sequence in file.get_sequence_iterator():
        for character_index, character in enumerate(sequence):
            if character_index >= read_length:
                break
            sequence_counts[character][character_index] += 1

    if not include_N:
        del sequence_counts["N"]

    return sequence_counts


def get_template_distances(template, FASTQ_file_path):
    """
    Given a FASTQ file path, return its distribution of distances to the given
    template.

    :param template: A template sequence
    :param FASTQ_file_path: Path to a FASTQ file
    :return: A list of distances
    """

    distances = []

    file = FASTQ_File(FASTQ_file_path)

    for sequence in file.get_sequence_iterator():

        distance = DNA_utils.get_template_distance(template, sequence)

        distances.append(distance)

    return distances


def get_matching_sequence_counts(
        extract_FASTQ_file_path,
        candidate_FASTQ_file_path=None,
        template=None,
        distance_threshold=None,
        quality_threshold=30,
        exclude_match=False,
        extract_indices=None):
    """
    Given a FASTQ file with sequences to match and a corresponding FASTQ file
    with sequences to extract, extract all the sequences where the matching
    sequence matches the given template.

    :param extract_FASTQ_file_path: Path to a FASTQ file of
        sequences to extract
    :param candidate_FASTQ_file_path: Path to a FASTQ file
    :param template: The template sequence that reads in the matching file
        should match
    :param distance_threshold: How far off from the template a transcript can
        be to still be considered matching.
    :param quality_threshold: The quality that all reads in the barcode UMI
        sequence must meet to be considered
    :param exclude_match: Whether to include (False) or exclude (True) the
        matches
    :param extract_indices: The indices to extract. If None, extracts all
    :return: A list of tuples of extract sequences and their count
    """

    line_index = 0
    extract_file = FASTQ_File(extract_FASTQ_file_path)
    extract_iterator = iter(extract_file.get_sequence_quality_iterator())

    if candidate_FASTQ_file_path is not None:
        candidate_file = FASTQ_File(candidate_FASTQ_file_path)
        candidate_iterator = iter(candidate_file.get_sequence_quality_iterator())
    else:
        candidate_file = None
        candidate_iterator = None

    sequence_counts = {}

    candidate = None

    while True:

        try:
            extract, extract_quality = next(extract_iterator)
            if extract_indices is not None:
                extract = extract[extract_indices[0]:extract_indices[1]]
        except StopIteration:
            break

        if candidate_iterator is not None:
            try:
                candidate, candidate_quality = next(candidate_iterator)
            except StopIteration:
                break

        if candidate is not None:
            distance = utils.get_sequence_distance(template, candidate)

            if distance > distance_threshold:
                if not exclude_match:
                    line_index += 1
                    continue
            elif exclude_match:
                line_index += 1
                continue

        meets_quality_threshold = True

        if quality_threshold is not None:

            quality_scores = \
                FASTQ_utils.convert_quality_string_to_quality_score(
                    extract_quality
                )

            if min(quality_scores) < quality_threshold:
                meets_quality_threshold = False

        if not meets_quality_threshold:
            line_index += 1
            continue

        if extract not in sequence_counts:
            sequence_counts[extract] = 1
        else:
            sequence_counts[extract] += 1

        line_index += 1

    if candidate_file is not None:
        candidate_file.close()
    extract_file.close()

    sequence_counts = [(sequence, count) for
                       sequence, count in sequence_counts.items()]

    return sequence_counts


def get_read_count(FASTQ_file_path):
    """
    Get the number of reads in a FASTQ file path. Line count / 4

    :param FASTQ_file_path: Path to a FASTQ file
    :return: The number of reads
    """

    file = FASTQ_File(FASTQ_file_path)

    read_count = 0

    for _ in file.get_quality_iterator():
        read_count += 1

    return read_count
