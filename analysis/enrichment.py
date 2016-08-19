import matplotlib.pyplot as plt
import seaborn as sns
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

def plot_distribution(enrichment_values):

    sns.distplot(enrichment_values,bins=100,kde=False,rug=False)
    plt.xlabel('log10(enrichment)')
    plt.ylabel('sequence counts')
    plt.show()

sns.set(color_codes=True)