import gzip


class FASTQ_File:

    def __init__(self, file_path):

        if file_path.endswith(".gz"):
            self._is_compressed = True
        else:
            self._is_compressed = False

        self._file_path = file_path
        self._file = None
        self._is_open = False

    def get_sequence_iterator(self):

        return FASTQ_Sequence_Iterator(self)

    def get_quality_iterator(self):

        return FASTQ_Quality_Iterator(self)

    def get_sequence_quality_iterator(self):

        return FASTQ_Sequence_Quality_Iterator(self)

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


class FASTQ_Sequence_Iterator:

    def __init__(self, FASTQ_file):
        self._file = FASTQ_file
        self._line_index = None

    def __iter__(self):
        self._file.open()
        self._next_line_index = 0
        return self

    def __next__(self):

        line = self._file.file_handle.readline()
        self._next_line_index += 1

        while self._next_line_index % 4 != 2:
            line = self._file.file_handle.readline()
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

        line = self._file.file_handle.readline()
        self._next_line_index += 1

        while self._next_line_index % 4 != 0:
            line = self._file.file_handle.readline()
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

        sequence_line = self._file.file_handle.readline()
        self._next_line_index += 1

        while self._next_line_index % 4 != 2:
            sequence_line = self._file.file_handle.readline()
            self._next_line_index += 1

        while self._next_line_index % 4 != 0:
            quality_line = self._file.file_handle.readline()
            self._next_line_index += 1

        if not sequence_line or not quality_line:
            self._file.close()
            raise StopIteration

        return sequence_line.strip(), quality_line.strip()

