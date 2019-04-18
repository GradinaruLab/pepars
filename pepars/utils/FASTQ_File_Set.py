from . import FASTQ_File


class FASTQ_File_Set:

    def __init__(self, file_paths):

        self._files = []

        for file_path in file_paths:
            self._files.append(FASTQ_File(file_path))
        self._is_open = False

    def get_sequence_iterator(self):

        return FASTQ_Set_Sequence_Iterator(self)

    def open(self):

        if self._is_open:
            raise ValueError("Can't open already open FASTQ files")

        for file in self._files:
            file.open()

        self._is_open = True

    def close(self):

        if not self._is_open:
            return

        for file in self._files:
            file.close()

        self._is_open = False

    @property
    def file_paths(self):
        return self._file_paths

    @property
    def files(self):
        return self._files


class FASTQ_Set_Sequence_Iterator:

    def __init__(self, FASTQ_file_set):
        self._file_set = FASTQ_file_set
        self._sequence_iterators = None

    def __iter__(self):

        self._sequence_iterators = []

        for file in self._file_set.files:
            iterator = file.get_sequence_iterator()
            self._sequence_iterators.append(iterator)
            iter(iterator)

        return self

    def __next__(self):

        sequences = []

        is_at_least_one_sequence_valid = False

        for iterator in self._sequence_iterators:

            try:
                sequence = next(iterator)
                is_at_least_one_sequence_valid = True
            except StopIteration:
                sequence = None

            sequences.append(sequence)

        if not is_at_least_one_sequence_valid:
            self._file_set.close()
            raise StopIteration

        return sequences
