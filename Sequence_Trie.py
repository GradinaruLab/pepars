from . import DNA


class Sequence_Trie_Node(object):

    def __init__(self, num_sequence_elements):

        self._children = [None] * num_sequence_elements
        self._leaf_value = None

    def add(self, sequence, value=None):

        node = self
        for element in sequence:

            if node._children[element] is None:
                node._children[element] = Sequence_Trie_Node(len(self._children))

            node = node._children[element]

        node._leaf_value = value

    def find(self, sequence):

        node = self

        for element in sequence:

            child_node = node._children[element]

            if child_node is None:
                return False
            else:
                node = child_node

        return True

    def get_value(self, sequence):

        node = self

        for element in sequence:

            child_node = node._children[element]

            if child_node is None:
                return None
            else:
                node = child_node

        return node._leaf_value


class Sequence_Trie(object):

    def __init__(self, by_nucleotide=True):

        self._element_index_map = {}

        if by_nucleotide:

            nucleotide_index = 0
            for nucleotide in DNA.get_nucleotides():
                self._element_index_map[nucleotide] = nucleotide_index
                nucleotide_index += 1
        else:
            amino_acid_index = 0
            for amino_acid in DNA.get_amino_acids():
                self._element_index_map[amino_acid] = amino_acid_index
                amino_acid_index += 1

        self._root = Sequence_Trie_Node(len(self._element_index_map))

    def convert_sequence_to_index_array(self, sequence):

        index_array = [0]*len(sequence)

        for element_index, element in enumerate(sequence):
            index_array[element_index] = self._element_index_map[element]

        return index_array

    def add(self, sequence, value=None):

        self._root.add(self.convert_sequence_to_index_array(sequence), value=value)

    def find(self, sequence):

        return self._root.find(self.convert_sequence_to_index_array(sequence))

    def get_value(self, sequence):

        return self._root.get_value(self.convert_sequence_to_index_array(sequence))

