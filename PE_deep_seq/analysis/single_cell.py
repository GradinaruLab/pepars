import numpy
import pandas

from ..utils import Sequence_Trie
from ..utils import DNA as DNA_utils


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


def add_barcode_UMI_counts_to_gene_counts(gene_counts_file_path,
                                          barcode_UMI_counts,
                                          gene_name,
                                          UMI_count_threshold=0):
    """
    Given a set of barcodes and their UMI counts, adds them to the given gene
    counts matrix.

    :param gene_counts_file_path: Path to a CSV file containing genes and their
        count in each barcode
    :param barcode_UMI_counts: A dictionary of barcodes, each entry of which is
        a dictionary of UMIs to their counts
    :param gene_name: The name of the gene being added
    :param UMI_count_threshold: How many counts of a UMI must exist to be
        considered valid
    :return: The gene_counts DataFrame
    """

    gene_counts = pandas.read_csv(gene_counts_file_path, index_col=0)

    valid_barcodes = [barcode[0:-2] for barcode in gene_counts.columns]

    barcodes_trie = Sequence_Trie(DNA_utils.NUCLEOTIDE_INDEX_MAP)

    for barcode in valid_barcodes:
        barcodes_trie.add(barcode)

    gene_counts.drop(gene_name, inplace=True, errors="ignore")

    gene_count_dict = {}

    for barcode, UMI_counts in barcode_UMI_counts.items():

        if "N" in barcode:
            continue

        if not barcodes_trie.find(barcode):
            continue

        UMI_count = 0

        for UMI, count in UMI_counts.items():
            if count > UMI_count_threshold:
                UMI_count += 1

        gene_count_dict[barcode + "-1"] = UMI_count

    new_gene_counts = pandas.DataFrame.from_dict(gene_count_dict,
                                                 orient="index",
                                                 columns=[gene_name])

    gene_counts = gene_counts.append(new_gene_counts.transpose(), sort=False)
    gene_counts = gene_counts.fillna(0)
    gene_counts = gene_counts.astype(dtype=numpy.uint32)

    gene_counts.to_csv(gene_counts_file_path, sep=",", encoding="utf-8",
                       chunksize=1000)

    return gene_counts


def get_barcodes_from_gene_counts(gene_counts_file_path):
    """
    Given a gene counts matrix, return all the barcodes

    :param gene_counts_file_path: Path to a CSV file containing genes and their
        count in each barcode
    :return: A list of cell barcodes
    """

    gene_counts = pandas.read_csv(gene_counts_file_path, index_col=0)

    valid_barcodes = [barcode[0:-2] for barcode in gene_counts.columns]

    return valid_barcodes
