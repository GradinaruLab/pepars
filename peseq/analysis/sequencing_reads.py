from ..utils import utils


def get_read_length(FASTQ_file_path):
    """
    Get the read length from a FASTQ file. Assumes all reads are of the same
    length

    :param FASTQ_file_path: Path to an uncompressed FASTQ file
    :return: The length of reads, as determined by the first read in the file
    """

    file = open(FASTQ_file_path)

    # Read the first two lines to get the length of the sequence
    file.readline()
    line = file.readline()

    read_length = len(line.strip())

    file.close()

    return read_length


def get_nucleotide_distribution(FASTQ_file_path):
    """
    Given a FASTQ file path, return its distribution of nucleotides in each
    position

    :param FASTQ_file_path: Path to an uncompressed FASTQ file
    :return: A dictionary of arrays, one entry for each unique character in the
        sequencing reads. Each array is the length of the reads, and the entries
        are the counts of that character at that position.
    """

    read_length = get_read_length(FASTQ_file_path)

    # Initialize the array - nucleotides + 1 for N
    sequence_counts = {}

    file = open(FASTQ_file_path)

    line_index = 0

    while True:
        line = file.readline()
        if not line:
            break
        if line_index % 4 == 1:
            for character_index, character in enumerate(line.strip()):

                if character not in sequence_counts:
                    sequence_counts[character] = [0] * read_length

                sequence_counts[character][character_index] += 1
        line_index += 1
    file.close()

    return sequence_counts


def get_template_distances(template, FASTQ_file_path):
    """
    Given a FASTQ file path, return its distribution of distances to the given
    template.

    :param template: A template sequence
    :param FASTQ_file_path: Path to an uncompressed FASTQ file
    :return: A list of distances
    """

    distances = []

    line_index = 0
    file = open(FASTQ_file_path)

    while True:
        line = file.readline()
        if not line:
            break
        if line_index % 4 == 1:
            sequence = line.strip()

            distance = utils.get_sequence_distance(template, sequence)

            distances.append(distance)
        line_index += 1
    file.close()

    return distances


def get_matching_sequence_counts(
        extract_FASTQ_file_path,
        candidate_FASTQ_file_path=None,
        template=None,
        distance_threshold=None,
        quality_threshold=30,
        exclude_match=False):
    """
    Given a FASTQ file with sequences to match and a corresponding FASTQ file
    with sequences to extract, extract all the sequences where the matching
    sequence matches the given template.

    :param extract_FASTQ_file_path: Path to an uncompressed FASTQ file of
        sequences to extract
    :param candidate_FASTQ_file_path: Path to an uncompressed FASTQ file
    :param template: The template sequence that reads in the matching file
        should match
    :param distance_threshold: How far off from the template a transcript can
        be to still be considered matching.
    :param quality_threshold: The quality that all reads in the barcode UMI
        sequence must meet to be considered
    :param exclude_match: Whether to include (False) or exclude (True) the
        matches
    :return: A dictionary of barcodes, each entry containing a dictionary of
        UMIs and the number of times this barcode/UMI combo appeared
    """

    line_index = 0
    extract_file = open(extract_FASTQ_file_path)
    if candidate_FASTQ_file_path is not None:
        candidate_file = open(candidate_FASTQ_file_path)
    else:
        candidate_file = None

    sequence_counts = {}

    candidate = None
    extract = None

    while True:

        extract_line = extract_file.readline()
        if not extract_line:
            break

        if candidate_file is not None:
            candidate_line = candidate_file.readline()
        else:
            candidate_line = None

        if line_index % 4 == 1:
            if candidate_file is not None:
                candidate = candidate_line.strip()
            extract = extract_line.strip()

        elif line_index % 4 == 3:

            quality_score_string = extract_line.strip()

            if candidate_file is not None:
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
                for char in quality_score_string:
                    if ord(char) - 33 < quality_threshold:
                        meets_quality_threshold = False
                        break

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

    :param FASTQ_file_path: Path to the uncompressed FASTQ file
    :return: The number of reads
    """

    file = open(FASTQ_file_path)

    line_count = 0

    while True:
        line = file.readline()
        if not line:
            break
        line_count += 1

    file.close()

    return int(line_count / 4)
