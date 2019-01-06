
def get_barcode_UMI_counts(sequence_counts):
    """
    Given a list of sequences and their counts, generate a dictionary of cell
    barcodes and their respective UMI counts. Assumes a 16:10 barcode/UMI split

    :param sequence_counts: A list of tuples of sequences and their count
    :return: A dictionary of sequences representing cell barcodes, each of which
        contains a dictionary of their UMIs and counts for each UMI
    """

    barcode_UMI_counts = {}

    for sequence, count in sequence_counts:

        barcode = sequence[0:16]
        UMI = sequence[16:]

        if barcode not in barcode_UMI_counts:
            barcode_UMI_counts[barcode] = {UMI: count}
        else:
            barcode_UMI_counts[barcode][UMI] = count

    return barcode_UMI_counts
