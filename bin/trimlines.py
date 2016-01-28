import sys

input_file = open(sys.argv[1], 'r')
output_file = open(sys.argv[2], 'w')

line_counter = 0

for line in input_file:
	if line_counter % 4 == 1:
		output_file.write('>Sequence' + str(int((line_counter + 3)/4)) + '\n')
		output_file.write(line)
	line_counter += 1

input_file.close()
output_file.close()