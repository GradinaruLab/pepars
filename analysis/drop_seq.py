import pandas
import h5py


def get_barcode_from_cell_name(cell_name):

    last_underscore = cell_name.rfind("_")
    last_hyphen = cell_name.rfind("-")

    return cell_name[last_underscore+1:last_hyphen]


def get_barcodes_from_transcript_matrix(file_path):

    barcodes = []

    if file_path[-4:] == ".csv":
        data_frame = pandas.read_csv(file_path)

        columns = list(data_frame.columns.values)

        barcode_counts = []

        # -1 to ignore the first column (since it's a header)
        num_cells = len(columns) - 1

        for cell_index in range(0, num_cells):

            cell_name = columns[cell_index + 1]

            barcode = get_barcode_from_cell_name(cell_name)

            barcodes.append(barcode)
            barcode_counts.append(0)
    elif file_path[-3:] == ".h5":
        file = h5py.File(file_path)
        transcriptome_name = list(file.keys())[0]
        barcodes = \
            [get_barcode_from_cell_name(x.decode("UTF-8"))
             for x in file[transcriptome_name]['barcodes']]

    return barcodes