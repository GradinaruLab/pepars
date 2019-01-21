import numpy
import math


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


def convert_sequence_label_dict_and_weights_to_matrices_and_weights(sequence_label_dict, sequence_weights):

    num_sequences = len(sequence_label_dict)
    num_sequence_elements = len(list(sequence_label_dict.items())[0][0])

    sequence_matrix = numpy.empty([num_sequences, num_sequence_elements], dtype=numpy.str)
    label_matrix = numpy.zeros([num_sequences], dtype=numpy.float32)
    weight_matrix = numpy.zeros([num_sequences, 1], dtype=numpy.float32)

    sequence_index = 0
    for sequence, feature in sequence_label_dict.items():
        for sequence_element_index in range(0, num_sequence_elements):
            sequence_matrix[sequence_index][sequence_element_index] = sequence[sequence_element_index]
        label_matrix[sequence_index] = feature
        weight_matrix[sequence_index] = sequence_weights[sequence]
        sequence_index += 1

    weight_matrix = weight_matrix.squeeze()
    return (sequence_matrix, label_matrix, weight_matrix)


def convert_sequence_label_dict_and_weights_to_matrices(sequence_label_dict, sequence_weights):

    num_sequences = len(sequence_label_dict)
    num_sequence_elements = len(list(sequence_label_dict.items())[0][0])

    sequence_matrix = []#numpy.empty([num_sequences, num_sequence_elements], dtype=numpy.str)
    label_matrix = []#numpy.zeros([num_sequences, 1], dtype=numpy.float32)

    for sequence, feature in sequence_label_dict.items():

        sequence_count_log = math.ceil(sequence_weights[sequence])

        for i in range(sequence_count_log):
            sequence_array = []
            for sequence_element_index in range(0, num_sequence_elements):
                sequence_array.append(sequence[sequence_element_index])

            sequence_matrix.append(sequence_array)
            label_matrix.append(feature)

    sequence_matrix = numpy.array(sequence_matrix, dtype=numpy.str)

    rng_state = numpy.random.get_state()
    numpy.random.shuffle(sequence_matrix)
    label_matrix = numpy.array(label_matrix, dtype=numpy.float32)
    label_matrix = numpy.reshape(label_matrix, (label_matrix.shape[0], 1))
    numpy.random.set_state(rng_state)
    numpy.random.shuffle(label_matrix)

    return (sequence_matrix, label_matrix)


def convert_string_keys_to_ints(dictionary):

    new_dictionary = {}

    for key, value in dictionary.items():
        new_dictionary[int(key)] = value

    return new_dictionary


def clean_string_for_filename(dirty_string):

    clean_string = dirty_string
    forbidden_characters = ["\\", "/", ":", "*", "?", "\"", "<", ">", "|"]
    for character in forbidden_characters:
        clean_string = clean_string.replace(character, "")
    return clean_string


def get_sequence_distance(sequence_1, sequence_2):

    min_length = min(len(sequence_1), len(sequence_2))
    max_length = max(len(sequence_1), len(sequence_2))

    distance = 0

    for character_index in range(0, min_length):
        if sequence_1[character_index] != sequence_2[character_index]:
            distance += 1

    distance += max_length - min_length

    return distance
