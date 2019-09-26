from ..utils.AminoAcid import AminoAcid
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import seaborn as sns
from ..utils import DNA
import pandas

# Usage
# enrichment_sequences_array: array of sequences that were under or above user's chosen enrichment threshold
# matrix property: amino acid property of interest. This function calculates the property for each amino acid in the sequence.
# output: creates an M x N matrix where M is the number of sequences and N is the number of positions in the sequence. Each value in the matrix represents corresponding property of interest value.
def generate_matrix_of_interest(enrichment_sequences_array,matrix_property):

    num_positions = len(enrichment_sequences_array[0])
    matrix_of_interest = np.zeros((len(enrichment_sequences_array),num_positions))
    for idx in range(0,len(enrichment_sequences_array)):
        current_sequence = enrichment_sequences_array[idx]
        for aa_idx in range(0,len(current_sequence)):
            current_amino_acid = AminoAcid(current_sequence[aa_idx])
            if (current_amino_acid.is_valid):
                    matrix_of_interest[idx][aa_idx]=current_amino_acid.properties[matrix_property]
            else:
                print(current_sequence)
    return matrix_of_interest


# Usage
# enrichment_sequences_array: array of sequences that were under or above user's chosen enrichment threshold
# matrix property: amino acid property of interest. This function calculates the property for the whole sequence.
# output: an M x 1 array where M is the number of sequences. The value at each position in the array corresponds to the property of interest.
def generate_array_of_interest(enrichment_sequences_array,array_property):
    array_of_interest = []
    for idx in range(0,len(enrichment_sequences_array)):
        current_sequence = enrichment_sequences_array[idx]
        if (str(array_property)=='aromaticity'):
            array_of_interest.append(ProteinAnalysis(str(current_sequence)).aromaticity())
        elif (str(array_property) == 'flexibility'):
            array_of_interest.append(ProteinAnalysis(str(current_sequence)).flexibility())
    return array_of_interest
    
def plot_amino_acid_property_distribution_from_matrix(matrix_of_interest,property_name_string,plot_title):
    
    global next_figure_index

    plt.figure(next_figure_index)
    #print next_figure_index
    next_figure_index += 1

    sequence_length = matrix_of_interest.shape[1]
    max_of_matrix = matrix_of_interest.max()
    min_of_matrix = matrix_of_interest.min()
    plt.title(plot_title)
    for i,v in enumerate(range(sequence_length)):
        v = v + 1
        ax1 = plt.subplot(sequence_length,1,v)
        ax1.set_xlim(min_of_matrix,max_of_matrix)
        sns.distplot(matrix_of_interest[:,i],bins=100,kde=False,rug=False,ax = ax1)
    plt.subplots_adjust(hspace=0.5)
    plt.xlabel(property_name_string)
    plt.title(plot_title)

def plot_amino_acid_property_distribution_from_array(array_of_interest_high,array_of_interest_low,property_name_string,plot_title_high,plot_title_low):
    # Determine minimum and maximum for the axis 
    max_of_array_high = max(array_of_interest_high)
    max_of_array_low =  max(array_of_interest_low)
    min_of_array_high = min(array_of_interest_high)
    min_of_array_low = min(array_of_interest_low)
    max_array = max(max_of_array_high,max_of_array_low)
    min_array = min(min_of_array_high,min_of_array_low)
    f,(ax1,ax2) = plt.subplots(2)
    sns.distplot(array_of_interest_high,bins=100,kde=False,ax=ax1)
    plt.xlabel(property_name_string)
    plt.ylabel('sequence counts')
    plt.title(plot_title_high)
    sns.distplot(array_of_interest_low,bins=100,kde=False,ax=ax2)
    plt.xlabel(property_name_string)
    plt.ylabel('sequence counts')
    plt.title(plot_title_low)
    plt.show(block=False)


def get_amino_acid_counts_by_position(sequences, counts=None):

    if isinstance(sequences, dict):
        counts = sequences.copy()
        sequences = list(sequences.keys())
        counts = [counts[x] for x in sequences]

    unique_AAs = DNA.get_amino_acids()
    unique_AAs = sorted(list(unique_AAs))

    counts_by_position = np.zeros((len(unique_AAs), len(sequences[0])))

    for sequence_index, sequence in enumerate(sequences):
        for character_index, character in enumerate(sequence):
            if counts is None:
                counts_by_position[unique_AAs.index(character), character_index] += 1
            else:
                counts_by_position[unique_AAs.index(character), character_index] += counts[sequence_index]

    counts_by_position = pandas.DataFrame(counts_by_position)
    counts_by_position.index = unique_AAs
    counts_by_position.columns = list(range(1, len(sequences[0]) + 1))
    return counts_by_position


def get_amino_acid_codon_biases(templates, template_ratios=None,
                                allow_stop_codons=False):

    if isinstance(templates, str):
        templates = [templates]

    template_length = len(templates[0])

    for template in templates:
        if len(template) != template_length:
            raise ValueError("All templates must have same length!")

    if template_length % 3 != 0:
        raise ValueError("Template must be a multiple of 3!")

    num_templates = len(templates)

    if template_ratios is None:
        template_ratios = np.ones((num_templates,)) / num_templates

    template_length = int(template_length / 3)

    unique_AAs = DNA.get_amino_acids()
    unique_AAs = sorted(list(unique_AAs))
    num_AAs = len(unique_AAs)

    weighted_AA_ratios = np.zeros((num_AAs, template_length))

    for template_index, template in enumerate(templates):

        amino_acid_ratios = np.zeros((num_AAs, template_length))

        for codon_index in range(template_length):

            for codon, amino_acid in DNA.CODON_AA_MAP.items():

                if not allow_stop_codons and amino_acid == "#":
                    continue

                is_codon_allowed = True
                for nucleotide_index in range(0, 3):
                    template_nucleotide_index = codon_index * 3 + nucleotide_index
                    allowable_nucleotides = DNA.IUPAC_GRAMMAR_MAP[template[template_nucleotide_index]]
                    if codon[nucleotide_index] not in allowable_nucleotides:
                        is_codon_allowed = False
                        break
                if is_codon_allowed:
                    amino_acid_ratios[unique_AAs.index(amino_acid), codon_index] += 1

        amino_acid_ratios /= amino_acid_ratios.sum(axis=0)
        weighted_AA_ratios += amino_acid_ratios * template_ratios[template_index]

    weighted_AA_ratios = pandas.DataFrame(weighted_AA_ratios, columns=list(range(1, template_length + 1)))
    weighted_AA_ratios.index = unique_AAs
    return weighted_AA_ratios

sns.set(color_codes=True)
next_figure_index = 1