from utils.AminoAcid import AminoAcid
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import seaborn as sns

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
                if (str(matrix_property) == 'molecular weight'):
                    matrix_of_interest[idx][aa_idx]=current_amino_acid.molecular_weight
                elif (str(matrix_property) == 'gravy'):
                    matrix_of_interest[idx][aa_idx]=current_amino_acid.hydrophobicity
            else:
                print current_sequence
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
    for i,v in enumerate(xrange(sequence_length)):
        v = v + 1
        ax1 = subplot(sequence_length,1,v)
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

sns.set(color_codes=True)
next_figure_index = 1