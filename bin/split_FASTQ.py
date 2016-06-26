# Reads in a FASTQ file and breaks it out into two corresponding files: one just
# containing the sequence, the other just containing the quality score
# Assigns a unique ID to each sequence

import sys

input_file = open(sys.argv[1], 'r')
fasta_file = open(sys.argv[2], 'w')
quality_file = open(sys.argv[3], 'w')

line_counter = 0
sequence_count = 0

for line in input_file:
	if line_counter % 4 == 1:
		sequence_count += 1
		fasta_file.write('>Sequence' + str(sequence_count) + '\n')
		fasta_file.write(line)
	elif line_counter % 4 == 3:
		quality_file.write('>Sequence' + str(sequence_count) + '\n')
		quality_file.write(line)
	line_counter += 1

input_file.close()
fasta_file.close()
quality_file.close()