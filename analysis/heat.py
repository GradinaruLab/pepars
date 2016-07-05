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


#Returns the number of combinations of 'acid' that endt with 'g' or 't'
def get_amino_factor(acid, gencode):
	count=0
	for key, value in gencode.iteritems():
		if(value.lower()==acid.lower() and (key[2].lower()=='t' or key[2].lower()=='g')):
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

class heatmap:
	def __init__(self, data = None, x_labels = None, y_labels = None, title = None):
		self.data = np.array(data)
		self.x_labels = x_labels
		self.y_labels = y_labels
		self.title = title
		self.start = None
 		self.midpoint = None
 		self.stop = None

	#Plots data according to 
	def normalized_sequence_counts(self, file, by_amino_acid = True, count_threshold = 10, filter_invalid = True):

#		get sequencedata
		slib = Sequence_Library(file);
		slib_count = slib.get_sequence_counts(by_amino_acid=by_amino_acid, count_threshold=count_threshold, filter_invalid=filter_invalid)

		sequence_matrix = []
		sequence_counts = []
		for sequence, count in slib_count.iteritems():
			sequence_matrix.append(list(sequence))
			sequence_counts.append(count)
		
		sequence_matrix = np.array(sequence_matrix)

#		hash different values of amino acids
#		SORT!
		acid_labels = {}
		acid_counts = {}
		for i, acid in enumerate(sorted(set(gencode.values()))):
			acid_labels[acid] = i
			acid_counts[acid] = get_amino_factor(acid, gencode)

		acids, positions = sequence_matrix.shape

#		biased_counts holds the counts for each amino acid for each position before normalization
		biased_counts = []

		for position in range(positions):
#			aminos is a column and represents the amino acids in position 'position'
			aminos = sequence_matrix[:,position] #column-wise
			row = [0 for acid in acid_labels]

#			note that aminos and sequence_counts have the same indices
			for sequence, amino in enumerate(aminos):
				if amino in acid_labels:
					row[acid_labels[amino]] += sequence_counts[sequence]
#					else: print 'acid: '+'\''+amino+'\''
			biased_counts.append(row)

#		Prepare the weight matrix
		weights = np.eye(len(acid_labels))
		for acid, count in acid_counts.iteritems():
#			get position of current acid in our enumeration			
			i = acid_labels[acid]
			weights[i, i] = 1.0/count
		
		'''print acid_labels
		print acid_counts
		for w in weights:
			print w
		'''
		self.x_labels = ['' for i in acid_labels]
		for acid, index in acid_labels.iteritems():
			self.x_labels[index]=acid

		self.y_labels = range(1,1+positions)
		self.data = np.dot(biased_counts, weights)/sum(sequence_counts)
		self.midpoint = 1.0/(32.0*np.amax(self.data)) if np.amax(self.data)>1.0/32.0 else 1.0
		self.start=np.amin(self.data)
		self.stop=1.0
		max=np.amax(self.data)
		if max<1.0/32.0:
			self.midpoint = 1.0
			self.stop=max*32.0

		return self
		'''
		self.midpoint = 1.0*np.sum(sequence_counts)/(32.0*(np.amax(normalized_data)))
		'''
		#heatmap(dataset, x_axis = x_axes, y_axis= y_axes, title = titles, midpoints = midpoints)
		#return dataset

	#Pre:	heatmap_objects is a heatmap object or a list of such objects
	#		Show is either Boolean or None.

	#Post:	The function has produced a heatmap for each 2D heatmap.
	#		If show is True, a heatmap has been drawn for each heatmap object, otherwise not

	@staticmethod
	def draw(heatmap_objects, show=True):
		if not type(heatmap_objects).__name__=='list':
			heatmap_objects=[heatmap_objects]
		for heatmap_object in heatmap_objects:
			assert(heatmap_object.__class__.__name__=='heatmap')
			data = heatmap_object.data
			row_labels = heatmap_object.y_labels if heatmap_object.y_labels else range(data.shape[0])
			column_labels = heatmap_object.x_labels if heatmap_object.x_labels else range(data.shape[1])
			title = heatmap_object.title if heatmap_object else ''
			midpoint = heatmap_object.midpoint if heatmap_object.midpoint else 0.5

			# set appropriate font and dpi
			fig, ax = plt.subplots()
			sns.set(font_scale=1.2)
			sns.set_style({"savefig.dpi": 100})

			
			start = heatmap_object.start if 0.0<=heatmap_object.start<=1.0 else 0.0
			stop =  heatmap_object.stop if 0.0<=heatmap_object.stop<=1.0 else 1.0

			ax = sns.heatmap(data,
				cmap=heatmap.shiftedColorMap(plt.cm.coolwarm,
					start=start,
					midpoint = midpoint,
					stop = stop),
				linewidths=.1)
			fig = ax.get_figure()


			# put the major ticks at the middle of each cell
			ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor = False)
			ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor = False)

			if(row_labels):
				ax.set_yticklabels(row_labels, minor = False,rotation=0)
			if(column_labels):
				ax.set_xticklabels(column_labels, minor = False)

			if(title):
				ax.set_title(title)
				fig.canvas.set_window_title('Heatmap for ' + title)

		if show:
			plt.show()
		
	@staticmethod
	def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
		'''
		returns a new colormap with its neutral value(center) at midpoint
		0.0<=start<midpoint<stop<=1.0
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
		if midpoint < stop:
			shift_index = np.hstack([
				np.linspace(0.0, midpoint, 128, endpoint=False), 
				np.linspace(midpoint, 1.0, 129, endpoint=True)
			])
		elif midpoint >= stop:
			shift_index = np.linspace(0.0, 1.0, 128, endpoint=True)
		for reg, shift in zip(reg_index, shift_index):
			r, g, b, a = cmap(reg)

			colordict['red'].append((shift, r, r))
			colordict['green'].append((shift, g, g))
			colordict['blue'].append((shift, b, b))
			colordict['alpha'].append((shift, a, a))

		newcmap = matplotlib.colors.LinearSegmentedColormap(name, colordict)
		plt.register_cmap(cmap=newcmap)

		return newcmap
	@staticmethod
	def compare_maps(map_a, map_b, title=None, x_labels=None, y_labels=None):
		assert(map_a.__class__.__name__=='heatmap')
		assert(map_b.__class__.__name__=='heatmap')
		assert(map_a.data.shape==map_b.data.shape)

		if not title:
			title = map_a.title + ' vs '+ map_b.title
		if not x_labels:
			x = map_a.x_labels
		if not y_labels:
			y = map_a.y_labels
		data = map_a.data-map_b.data

		return heatmap(data=data, x_labels=x, y_labels=y, title=title)