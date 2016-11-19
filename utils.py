import numpy

def convert_string_to_char_array(string):
	char_array = []

	if type(string)!=str:
		for char in string:
			char_array.append(char);
		return char_array
	return list(string)

def convert_sequence_label_dict_to_matrices(sequence_label_dict):

	num_sequences = len(sequence_label_dict)
	num_sequence_elements = len(list(sequence_label_dict.items())[0][0])

	sequence_matrix = numpy.empty([num_sequences, num_sequence_elements], dtype=numpy.str)
	label_matrix = numpy.zeros([num_sequences, 1], dtype=numpy.float32)

	sequence_index = 0
	for sequence, feature in sequence_label_dict.items():
		for sequence_element_index in range(0, num_sequence_elements):
			sequence_matrix[sequence_index][sequence_element_index] = sequence[sequence_element_index]
			label_matrix[sequence_index] = feature
		sequence_index += 1

	return (sequence_matrix, label_matrix)

def convert_string_keys_to_ints(dictionary):

	new_dictionary = {}

	for key, value in dictionary.items():
		new_dictionary[int(key)] = value

	return new_dictionary
