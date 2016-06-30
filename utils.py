

def convert_string_to_char_array(string):
	char_array = []

	if type(string)!=str:
		for char in string:
			char_array.append(char);
		return char_array
	return list(string)

def convert_sequence_label_dict_to_matrices(sequence_label_dict):

	sequence_matrix = []
	label_matrix = []

	for sequence, feature in sequence_label_dict.items():
		sequence_matrix.append(convert_string_to_char_array(sequence))
		label_matrix.append([feature])

	return (sequence_matrix, label_matrix)