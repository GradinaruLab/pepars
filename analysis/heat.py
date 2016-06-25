import glob
import sys
import numpy as np
import os
from utils.Analysis_Set import Analysis_Set
import matplotlib
import matplotlib.pyplot as plt

from utils.Sequence_Library import Sequence_Library
from utils import utils
from utils import DNA as DNA
from analysis import ML_analysis
from fileio import sequence_data_parsing
import csv

#Returns the number of combinations of 'acid' that don't end 
def get_factor(acid, gencode):
	count=0
	for key, value in gencode.iteritems():
		if(value==acid and (key[2].lower()=='t' or key[2].lower()=='g')):
			count+=1
	return count

gencode = {
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
		'TAC':'Y', 'TAT':'Y', 'TAA':'#', 'TAG':'#', 'TGC':'C', 'TGT':'C', 'TGA':'#', 'TGG':'W'}

#Plots data according to 
def plot_data(files, by_amino_acid = True, count_threshold = 10, filter_invalid = True):
	if not type(files) == list:
		files = [files]
	dataset = []
	for sequence_data in files:

		slib = Sequence_Library(sequence_data);
		slib_count = slib.get_sequence_counts(by_amino_acid=by_amino_acid, count_threshold=count_threshold, filter_invalid=filter_invalid)

		sequence_matrix = []
		sequence_counts = []
		for sequence, count in slib_count.iteritems():
			sequence_matrix.append(list(sequence))
			sequence_counts.append(count)

		sequence_matrix = np.array(sequence_matrix)

		#hash different values of amino acids
		acid_labels = {}
		acid_counts = {}
		for i, acid in enumerate(set(gencode.values())):
			acid_labels[acid] = i
			acid_counts[acid] = get_factor(acid, gencode)

		acids, positions = sequence_matrix.shape

		#biased_counts holds the counts for each amino acid for each position before normalization
		biased_counts = []

		for position in range(positions):
			#aminos is a column and represents the amino acids in position 'position'
			aminos = sequence_matrix[:,position]
			row = [0 for acid in acid_labels]
			#note that aminos and sequence_counts have the same indices
			for sequence, amino in enumerate(aminos):
				if amino in acid_labels:
					row[acid_labels[amino]] += sequence_counts[sequence]
				else: print 'acid: '+'\''+amino+'\''
			biased_counts.append(row)

		#Prepare the weight matrix
		weights = np.eye(len(acid_labels))
		
		for acid, count in acid_counts.iteritems():
			i = acid_labels[acid] #get position of current acid in our enumeration
			weights[i, i] = 1/count


		fig_labels_acids = ['' for i in acid_labels]
		for acid, index in acid_labels.iteritems():
			fig_labels_acids[index]=acid

		data = np.dot(biased_counts, weights)
		dataset.append(data)

		row_labels = fig_labels_acids
		column_labels = range(1,1+positions)

		fig, ax = plt.subplots()
		fig.title='Hell nah'
		heatmap = ax.pcolor(data)

		# put the major ticks at the middle of each cell
		ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor = False)
		ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor = False)

		ax.set_xticklabels(row_labels, minor = False)
		ax.set_yticklabels(column_labels, minor = False)

		ax.set_title(sequence_data)
		fig.canvas.set_window_title('Heatmap for ' + sequence_data) 
	plt.show()
	return dataset


'''
files = ['Data/21_perfect_match_sequence_counts.csv',
		'Data/22_perfect_match_sequence_counts.csv',
		'Data/23_perfect_match_sequence_counts.csv']
dataset =  plot_data(files)
'''