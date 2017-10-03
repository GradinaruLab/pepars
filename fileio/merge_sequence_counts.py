import sys
from collections import defaultdict
import csv

num_files = len(sys.argv) - 2

files = []

for file_index in range(0,num_files):
    files.append(csv.reader(open(sys.argv[file_index + 1], 'r',encoding='utf8',newline='')))

variant_sequences = defaultdict(list)

writer = csv.writer(open(sys.argv[num_files + 1], 'w',encoding='utf8',newline=''))

header_row = ['Sequence']

for file_index in range(0, num_files):

    header_row.append(sys.argv[file_index+1])

    file_reader = files[file_index]

    for sequence_row in file_reader:

        sequence = sequence_row[0]

        num_list_entries = len(variant_sequences[sequence])

        for fill_in_index in range(num_list_entries, file_index + 1):
            variant_sequences[sequence].append(0)

        variant_sequences[sequence][file_index] = int(sequence_row[1])

writer.writerow(header_row)

for sequence, sequence_counts in sorted(variant_sequences.items(), key=lambda kv: sum(kv[1]), reverse=True):

    num_list_entries = len(sequence_counts)

    if num_list_entries != num_files:
        num_missing_list_entries = num_files - num_list_entries
        for i in range(0,num_missing_list_entries):
            sequence_counts.append(0)

    row = [sequence]

    for sequence_count in sequence_counts:
        row.append(sequence_count)

    writer.writerow(row)