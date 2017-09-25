import pandas
import os

def get_barcodes_from_transcript_matrix(file_path):

    barcodes = []

    file_path = os.path.join("..", "workspace", "drop_seq_1", "H2-small_matrix.csv")

    data_frame = pandas.read_csv(file_path)

    columns = list(data_frame.columns.values)

    barcode_counts = []

    # -1 to ignore the first column (since it's a header)
    num_cells = len(columns) - 1

    for cell_index in range(0, num_cells):

        cell_name = columns[cell_index + 1]

        last_underscore = cell_name.rfind("_")
        last_hyphen = cell_name.rfind("-")

        barcode = cell_name[last_underscore+1:last_hyphen]

        barcodes.append(barcode)
        barcode_counts.append(0)

    print(data_frame.sum())


    # for index, row in data_frame.iterrows():

    #     print(index)

    #     if index == 0:
    #         continue

    #     for cell_index in range(0, num_cells):

    #         count = row[cell_index + 1]

    #         barcode_counts[cell_index] += count

    # for barcode_index in range(0, len(barcodes)):

    #     if barcode_counts[barcode_index] == 0:
    #         num_zero_counts += 1

    # print("Num zero counts: %i" % num_zero_counts)

    return barcodes