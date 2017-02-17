import seaborn
import numpy
import matplotlib.pyplot as plt

def plot_sequence_count_histogram(sequence_counts):

    seaborn.distplot(numpy.log10(sequence_counts))
    plt.show()