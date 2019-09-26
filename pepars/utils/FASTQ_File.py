import gzip
import os
from enum import Enum


class FASTQ_File:

    def __init__(self, file_path):

        if file_path.endswith(".gz"):
            self._is_compressed = True
        elif not os.path.isfile(file_path) and \
                os.path.isfile(file_path + ".gz"):

            file_path += ".gz"
            self._is_compressed = True
        else:
            self._is_compressed = False

        self._file_path = file_path
        self._file = None
        self._is_open = False
        self._file_name = os.path.basename(file_path)

    def get_sequence_iterator(self):

        return FASTQ_Sequence_Iterator(self)

    def get_quality_iterator(self):

        return FASTQ_Quality_Iterator(self)

    def get_sequence_quality_iterator(self):

        return FASTQ_Sequence_Quality_Iterator(self)

    def get_read_count(self):

        read_count = 0

        for _ in self.get_sequence_iterator():
            read_count += 1

        return read_count

    def open(self):

        if self._is_open:
            raise ValueError("Can't open already open FASTQ file")

        if self._is_compressed:
            self._file = gzip.open(self._file_path, mode="rt", encoding="utf-8")
        else:
            self._file = open(self._file_path, mode="r")

        self._is_open = True

    def close(self):

        if not self._is_open:
            return

        self._is_open = False

        self._file.close()

    @property
    def is_compressed(self):
        return self._is_compressed

    @property
    def file_path(self):
        return self._file_path

    @property
    def file_handle(self):
        return self._file

    @property
    def file_name(self):
        return self._file_name


class FASTQ_Sequence_Iterator:

    def __init__(self, FASTQ_file):
        self._file = FASTQ_file
        self._line_index = None

    def __iter__(self):
        self._file.open()
        self._next_line_index = 0
        return self

    def __next__(self):

        try:
            line = self._file.file_handle.readline()
        except OSError:
            self._file.close()
            raise StopIteration

        self._next_line_index += 1

        while self._next_line_index % 4 != 2:
            try:
                line = self._file.file_handle.readline()
            except OSError:
                self._file.close()
                raise StopIteration
            self._next_line_index += 1

        if not line:
            self._file.close()
            raise StopIteration

        return line.strip()


class FASTQ_Quality_Iterator:

    def __init__(self, FASTQ_file):
        self._file = FASTQ_file
        self._line_index = None

    def __iter__(self):
        self._file.open()
        self._next_line_index = 0
        return self

    def __next__(self):

        try:
            line = self._file.file_handle.readline()
        except OSError:
            self._file.close()
            raise StopIteration

        self._next_line_index += 1

        while self._next_line_index % 4 != 0:
            try:
                line = self._file.file_handle.readline()
            except OSError:
                self._file.close()
                raise StopIteration
            self._next_line_index += 1

        if not line:
            self._file.close()
            raise StopIteration

        return line.strip()


class FASTQ_Sequence_Quality_Iterator:

    def __init__(self, FASTQ_file):
        self._file = FASTQ_file
        self._line_index = None

    def __iter__(self):
        self._file.open()
        self._next_line_index = 0
        return self

    def __next__(self):

        try:
            sequence_line = self._file.file_handle.readline()
        except OSError:
            self._file.close()
            raise StopIteration

        self._next_line_index += 1

        while self._next_line_index % 4 != 2:
            try:
                sequence_line = self._file.file_handle.readline()
            except OSError:
                self._file.close()
                raise StopIteration
            self._next_line_index += 1

        while self._next_line_index % 4 != 0:
            try:
                quality_line = self._file.file_handle.readline()
            except OSError:
                self._file.close()
                raise StopIteration
            self._next_line_index += 1

        if not sequence_line or not quality_line:
            self._file.close()
            raise StopIteration

        return sequence_line.strip(), quality_line.strip()


class Illumina_FASTQ_File(FASTQ_File):

    def __init__(self, file_path):
        super().__init__(file_path)

        file_name = os.path.basename(file_path)

        file_name_parts = file_name.split("_")

        if len(file_name_parts) < 5:
            raise ValueError("Illumina FASTQ file should have at least 5 \
                             underscore-separated components!")

        self._sample_name = "_".join(file_name_parts[0:-4])
        self._sample_number = int(file_name_parts[-4][1:])
        self._lane_number = int(file_name_parts[-3][1:])

        read_type = file_name_parts[-2]

        if read_type[0] == "I":
            self._is_index_read = True
        else:
            self._is_index_read = False

        self._read_number = int(read_type[1:])
        self._file_number = int(file_name_parts[-1].split(".")[0])

    @property
    def sample_name(self):
        return self._sample_name

    @property
    def sample_number(self):
        return self._sample_number

    @property
    def lane_number(self):
        return self._lane_number

    @property
    def is_index_read(self):
        return self._is_index_read

    @property
    def read_number(self):
        return self._read_number

    @property
    def file_number(self):
        return self._file_number
