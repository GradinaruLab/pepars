import queue

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

    def __iter__(self):

        self._current_node_chain = [self._root]
        self._current_node_element_index = [0]

        return self

    def __next__(self):

        # Get the last node in the chain
        current_node = self._current_node_chain[-1]
        current_child_index = self._current_node_element_index[-1]

        next_child_index = None

        for child_index in range(current_child_index + 1,
                                 len(current_node.children)):
            if current_node.children[child_index] is not None:
                next_child_index = child_index + 1
                self._current_node_chain.append(
                    current_node.children[child_index])
                self._current_node_element_index.append(0)
            break

        # If we didn't find a next child,
        if next_child_index is None:
            if len(self._current_node_chain) == 1:
                raise StopIteration
            else:
                pass
