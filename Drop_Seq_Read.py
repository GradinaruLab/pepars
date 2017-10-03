
class Drop_Seq_Read(object):

    def __init__(self):
        self._cell_barcode = None
        self._UMI = None
        self._virus_sequence = None

    @property
    def cell_barcode(self):
        return self._cell_barcode

    @property
    def UMI(self):
        return self._UMI

    @property
    def virus_sequence(self):
        return self._virus_sequence

    @cell_barcode.setter
    def cell_barcode(self, new_cell_barcode):
        self._cell_barcode = new_cell_barcode

    @UMI.setter
    def UMI(self, new_UMI):
        self._UMI = new_UMI

    @virus_sequence.setter
    def virus_sequence(self, new_virus_sequence):
        self._virus_sequence = new_virus_sequence