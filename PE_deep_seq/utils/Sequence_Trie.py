from . import DNA as DNA_utils


class Sequence_Trie_Node(object):

    def __init__(self, element_index_map):

        self._children = [None] * len(element_index_map)
        self._leaf_value = None
        self._element_index_map = element_index_map

    def add(self, sequence, value=None):

        node = self
        for element in sequence:

            index = self._element_index_map[element]

            if node.children[index] is None:
                node.children[index] = \
                    Sequence_Trie_Node(self._element_index_map)

            node = node._children[index]

        node._leaf_value = value

    def find(self, sequence):

        node = self

        for element in sequence:

            child_node = node.children[self._element_index_map[element]]

            if child_node is None:
                return False
            else:
                node = child_node

        return True

    def get_value(self, sequence):

        node = self

        for element in sequence:

            child_node = node.children[self._element_index_map[element]]

            if child_node is None:
                return None
            else:
                node = child_node

        return node._leaf_value

    def get_child(self, sequence):

        node = self

        for element in sequence:

            child_node = node.children[self._element_index_map[element]]

            if child_node is None:
                return None
            else:
                node = child_node

        return node

    @property
    def children(self):
        return self._children


class Sequence_Trie(object):

    def __init__(self, by_nucleotide=True, allow_invalid=False):

        if by_nucleotide:
            self._element_index_map = dict(DNA_utils.NUCLEOTIDE_INDEX_MAP)
            if allow_invalid:
                self._element_index_map["N"] = \
                    len(DNA_utils.NUCLEOTIDE_INDEX_MAP)
        else:
            self._element_index_map = dict(DNA_utils.AMINO_ACID_INDEX_MAP)
            if allow_invalid:
                self._element_index_map["#"] = \
                    len(DNA_utils.AMINO_ACID_INDEX_MAP)

        self._root = Sequence_Trie_Node(self._element_index_map)

    def add(self, sequence, value=None):

        self._root.add(sequence, value=value)

    def find(self, sequence):

        return self._root.find(sequence)

    def get_value(self, sequence):

        return self._root.get_value(sequence)

    def get_node(self, sequence):

        return self._root.get_child(sequence)
