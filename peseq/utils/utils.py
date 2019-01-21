

def get_sequence_distance(sequence_1, sequence_2):

    min_length = min(len(sequence_1), len(sequence_2))
    max_length = max(len(sequence_1), len(sequence_2))

    distance = 0

    for character_index in range(0, min_length):
        if sequence_1[character_index] != sequence_2[character_index]:
            distance += 1

    distance += max_length - min_length

    return distance
