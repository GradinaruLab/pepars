import sys
import os
import time
from sklearn.metrics import roc_curve, auc
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import csv
import numpy
import random
import math
import time
from sklearn.neighbors.kde import KernelDensity

num_sample_range = []
num_sample_range.extend(list(range(1, 9, 1)))
# num_sample_range.extend(list(range(10, 90, 10)))
# num_sample_range.extend(list(range(100, 900, 100)))
# num_sample_range.extend(list(range(1000, 9000, 1000)))
# num_sample_range.extend(list(range(10000, 90000, 10000)))
# num_sample_range.extend(list(range(100000, 900000, 100000)))
#num_sample_range.extend(list(range(1000000, 9000000, 1000000)))
#num_sample_range.extend(list(range(10000000, 90000000, 10000000)))
sequence_length = 1
sequence_options = ['A', 'T', 'C', 'G']
num_bins = 10
weighted = True

plt.ion()
plt.show()
plt.figure(figsize=(15,15))

for num_samples in num_sample_range:

	sequence_counts = {}

	for sample_index in range(0, num_samples):

		sequence = ''

		for sequence_index in range(0, sequence_length):

			sequence += random.choice(sequence_options)


		if sequence not in sequence_counts:
			sequence_counts[sequence] = 1
		else:
			sequence_counts[sequence] += 1

	data_values = {}

	cdf_values = []

	for sequence, sequence_count in sequence_counts.items():

		if sequence_count not in data_values:
			if weighted:
				data_values[sequence_count] = sequence_count
				cdf_values.extend([sequence_count] * sequence_count)
			else:
				data_values[sequence_count] = 1
		else:
			if weighted:
				data_values[sequence_count] += sequence_count
				cdf_values.extend([sequence_count] * sequence_count)
			else:
				data_values[sequence_count] += 1

	total_value = sum(data_values.values())

	num_values = len(data_values)
	x_values = numpy.zeros([num_values], dtype=numpy.float32)
	y_values = numpy.zeros([num_values], dtype=numpy.float32)

	total_count = sum(data_values.values())
	print("Total count: %i" % total_count)
	for x, y in data_values.items():
		if x == 1:
			num_single_counts = y
			break

	print("Num single counts: %i" % num_single_counts)

	num_possible_counts = pow(len(sequence_options), sequence_length)

	print("Num possible counts: %i" % num_possible_counts)

	num_unique = len(sequence_counts)

	print("Num unique sampled: %i" % num_unique)

	num_unseen = num_possible_counts - num_unique

	print("Percent unseen: %.4f%%" % (num_unseen / num_possible_counts * 100))

	turing_number = num_single_counts / total_count
	print("Turing probability: %.4f%%" % (turing_number * 100))

	# index = 0
	# for x, y in data_values.items():
	# 	x_values[index] = x
	# 	y_values[index] = y / total_value
	# 	index += 1



	# plt.scatter(x_values, y_values, s=100)
	# plt.gca().set_xlim(left=0)
	# plt.title("Sequence counts from " + str(num_samples) + " random sequences of length " + str(sequence_length))
	# plt.draw()
	# plt.pause(2)
	# plt.clf()