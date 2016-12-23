import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from utils.AminoAcid import AminoAcid

def split_by_enrichment(sequence_enrichments, enrichment_threshold):

    num_positions = len(sequence_enrichments.keys()[0])

    above_enrichment_sequences = []
    below_enrichment_sequences = []
    enrichment_values=[]
    invalid_amino_acid = False;
    for sequence, enrichment in sequence_enrichments.items():
        for seq_idx in range(0,num_positions):
            if (AminoAcid(sequence[seq_idx]).is_valid):
                invalid_amino_acid = False
            else:
                invalid_amino_acid = True
                break
        if (invalid_amino_acid==False):
            enrichment_values.append(enrichment)
            if (enrichment > enrichment_threshold):
                above_enrichment_sequences.append(sequence)
            elif(enrichment < enrichment_threshold):
                below_enrichment_sequences.append(sequence)

    return above_enrichment_sequences, below_enrichment_sequences


# note: std_dev_above_factor can be negative if the user wants to find a threshold below the mean
def get_threshold(vector_of_interest,std_dev_above_factor=1):
    #Calculate the mean and stdev of the enrichment values
    vector_mean = np.mean(vector_of_interest)
    standard_dev = np.std(vector_of_interest)
    threshold = vector_mean + std_dev_above_factor*standard_dev
    return threshold


def plot_distribution(enrichment_values,threshold):

    sns.distplot(enrichment_values,bins=100,kde=False,rug=False)
    plt.xlabel('log10(enrichment)')
    plt.ylabel('sequence counts')
    plt.axvline(threshold)
    plt.show()

sns.set(color_codes=True)