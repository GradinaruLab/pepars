import glob
import sys
import numpy as np
import os
from utils.Analysis_Set import Analysis_Set
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

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
def normalized_heatmap(files, by_amino_acid = True, count_threshold = 10, filter_invalid = True):
	if not type(files) == list:
		files = [files]
	dataset = []
	x_axes = []
	y_axes = []
	titles = []
	midpoints = []
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

		
		normalized_data = np.dot(biased_counts, weights)

		x_axes.append(fig_labels_acids)
		y_axes.append(range(1,1+positions))
		titles.append(sequence_data)
		#midpoints.append(np.sum(sequence_counts))
		midpoint = 1.0*np.sum(sequence_counts)/(32.0*(np.amax(normalized_data)))
		midpoints.append(midpoint)
		dataset.append(normalized_data)

	heatmap(dataset, x_axis = x_axes, y_axis= y_axes, title = titles, midpoints = midpoints)
	return dataset

#Pre:	dataset is an array of one or more 2D arrays containing plottable data
#		x_axis is an array of one or more 1D arrays containing labels for the x-axis
#		y_axis is an array of one or more 1D arrays containing labels for the y-axis
#		titles is an array of one or more 1D arrays containing the name of each plot
#		no_seq is an array of one or more 1D arrays containing the number of sequences of each dataset
#		all parameters need to have the same length (i.e. number of heatmaps)

#Post:	The function has produced a heatmap for each 2D array in dataset
def heatmap(dataset, x_axis = None, y_axis = None, title = None, midpoints=None):
	number_of_heatmaps = len(dataset)
	if not x_axis:
		x_axis = [[] in range(number_of_heatmaps)]
	if not y_axis:
		y_axis = [[] in range(number_of_heatmaps)]
	if not title:
		title = [[] in range(number_of_heatmaps)]
	if not midpoints:
		midpoints = [[] in range(number_of_heatmaps)]
	for matrix, x, y, name, midpoint in zip(np.array(dataset), x_axis, y_axis, title, midpoints):
		data = matrix
		row_labels = x
		column_labels = y
		
		# set appropriate font and dpi
		fig, ax = plt.subplots()
		sns.set(font_scale=1.2)
		sns.set_style({"savefig.dpi": 100})
		
		ax = sns.heatmap(data, cmap=shiftedColorMap(plt.cm.coolwarm, midpoint = midpoint if 0.0<midpoint<1.0 else 0.5), linewidths=.1)
		fig = ax.get_figure()


		# put the major ticks at the middle of each cell
		ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor = False)
		ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor = False)

		if(row_labels):
			ax.set_xticklabels(row_labels, minor = False)
		if(column_labels):
			ax.set_yticklabels(column_labels, minor = False)

		if(name):
			ax.set_title(name)
			fig.canvas.set_window_title('Heatmap for ' + name)
		
	plt.show()
	

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
	'''
	returns a new colormap with its neutral value(center) at midpoint
	start=0.0<midpoint<stop=1.0
	the returned colormap is based on cmap
	'''
	colordict = {
		'red': [],
		'green': [],
		'blue': [],
		'alpha': []
	}

	# regular index to compute the colors
	reg_index = np.linspace(start, stop, 257)

	# shifted index to match the data
	shift_index = np.hstack([
		np.linspace(0.0, midpoint, 128, endpoint=False), 
		np.linspace(midpoint, 1.0, 129, endpoint=True)
	])

	for reg, shift in zip(reg_index, shift_index):
		r, g, b, a = cmap(reg)

		colordict['red'].append((shift, r, r))
		colordict['green'].append((shift, g, g))
		colordict['blue'].append((shift, b, b))
		colordict['alpha'].append((shift, a, a))

	newcmap = matplotlib.colors.LinearSegmentedColormap(name, colordict)
	plt.register_cmap(cmap=newcmap)

	return newcmap

def compare_maps(data_a, data_b, x_axis=None, y_axis=None, title = None, midpoints=None):
	data = np.array(data_a)-np.array(data_b)
	heatmap([data], x_axis = x_axis, y_axis = y_axis, title = title, midpoints=midpoints)