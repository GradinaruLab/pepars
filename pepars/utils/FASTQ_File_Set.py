import os
from . import FASTQ_File
from . import Illumina_FASTQ_File


class FASTQ_File_Set:

    def __init__(self, file_paths):

        self._files = []

        for file_path in file_paths:
            self._files.append(FASTQ_File(file_path))

        self._is_open = False

    def get_sequence_iterator(self):

        return FASTQ_Set_Sequence_Iterator(self)

    def get_sequence_quality_iterator(self):

        return FASTQ_Set_Sequence_Quality_Iterator(self)

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

    def get_read_count(self):

        read_count = self._files[0].get_read_count()

        return read_count

    @property
    def num_files(self):
        return len(self._files)

    @property
    def file_paths(self):
        return [file.file_path for file in self._files]

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

        self._is_open = True

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


class FASTQ_Set_Sequence_Quality_Iterator:

    def __init__(self, FASTQ_file_set):
        self._file_set = FASTQ_file_set
        self._sequence_quality_iterators = None

    def __iter__(self):

        self._sequence_quality_iterators = []

        for file in self._file_set.files:
            iterator = file.get_sequence_quality_iterator()
            self._sequence_quality_iterators.append(iterator)
            iter(iterator)

        self._is_open = True

        return self

    def __next__(self):

        sequences = []
        quality_scores = []

        is_at_least_one_sequence_valid = False

        for iterator in self._sequence_quality_iterators:

            try:
                sequence, quality_score = next(iterator)
                is_at_least_one_sequence_valid = True
            except StopIteration:
                sequence = None
                quality_score = None

            sequences.append(sequence)
            quality_scores.append(quality_score)

        if not is_at_least_one_sequence_valid:
            self._file_set.close()
            raise StopIteration

        return sequences, quality_scores


class Illumina_FASTQ_Set_Sequence_Iterator:

    def __init__(self, Illumina_FASTQ_file_set):

        self._file_set = Illumina_FASTQ_file_set
        self._sequence_iterators = None

    def __iter__(self):

        self._initiate_lane_file_number(0)

        self._is_open = True

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

            for file in self._file_set.files_by_lane_file_number[
                    self._lane_file_number]:
                file.close()

            next_lane_file_index = self._lane_file_index + 1

            if next_lane_file_index >= len(self._file_set.lane_file_numbers):
                raise StopIteration

            self._initiate_lane_file_number(next_lane_file_index)

            return self.__next__()

        return sequences

    def _initiate_lane_file_number(self, lane_file_index):

        self._lane_file_index = lane_file_index
        self._lane_file_number = self._file_set.lane_file_numbers[
            self._lane_file_index]

        self._sequence_iterators = []

        for file in self._file_set.files_by_lane_file_number[
                self._lane_file_number]:
            iterator = file.get_sequence_iterator()
            self._sequence_iterators.append(iterator)
            iter(iterator)


class Illumina_FASTQ_File_Set(FASTQ_File_Set):

    def __init__(self, FASTQ_directory_path, sample_name):
        super().__init__([])

        directory_entries = os.listdir(FASTQ_directory_path)

        # Files by unique sample/lane numbers, in case this sample is split
        # over multiple lanes
        self._files_by_lane_file_number = {}

        for entry in directory_entries:

            if entry.startswith(sample_name):

                file_path = os.path.join(FASTQ_directory_path, entry)
                FASTQ_file = Illumina_FASTQ_File(file_path)

                if FASTQ_file.sample_name != sample_name:
                    continue

                sample_lane_file_number = \
                    (FASTQ_file.sample_number,
                     FASTQ_file.lane_number,
                     FASTQ_file.file_number)

                if sample_lane_file_number not in self._files_by_lane_file_number:
                    self._files_by_lane_file_number[sample_lane_file_number] = []

                self._files_by_lane_file_number[sample_lane_file_number].\
                    append(FASTQ_file)

        for sample_lane_file_number, illumina_files in \
                self._files_by_lane_file_number.items():

            # Now we add the Illumina files in the following order: indexes,
            # followed by reads

            latest_read_number = 0

            self._files_by_lane_file_number[sample_lane_file_number] = []

            added_file = True

            while added_file:
                added_file = False
                for file in illumina_files:
                    if file.is_index_read and \
                            file.sample_number == sample_lane_file_number[0] and \
                            file.lane_number == sample_lane_file_number[1] and \
                            file.read_number == latest_read_number + 1:
                        self._files.append(file)
                        self._files_by_lane_file_number[sample_lane_file_number].\
                            append(file)
                        latest_read_number += 1
                        added_file = True
                        break

            latest_read_number = 0

            added_file = True

            while added_file:
                added_file = False
                for file in illumina_files:
                    if not file.is_index_read and \
                            file.sample_number == sample_lane_file_number[0] and \
                            file.lane_number == sample_lane_file_number[1] and \
                            file.read_number == latest_read_number + 1:
                        self._files.append(file)
                        self._files_by_lane_file_number[sample_lane_file_number].\
                            append(file)
                        latest_read_number += 1
                        added_file = True
                        break

        self._lane_file_numbers = \
            sorted(list(self._files_by_lane_file_number.keys()))

    def get_index_file(self, index_number=1):

        for file in self._files:
            if file.is_index_read and file.read_number == index_number:
                return file

        return None

    def get_read_file(self, read_number=1):

        for file in self._files:
            if not file.is_index_read and file.read_number == read_number:
                return file

        return None

    def get_sequence_iterator(self):

        return Illumina_FASTQ_Set_Sequence_Iterator(self)

    def get_sequence_quality_iterator(self):

        return FASTQ_Set_Sequence_Quality_Iterator(self)

    @property
    def lane_file_numbers(self):
        return self._lane_file_numbers

    @property
    def files_by_lane_file_number(self):
        return self._files_by_lane_file_number
