import csv
import pandas


def write_csv_file(file_name, header, rows):

    file = open(file_name, 'w')
    writer = csv.writer(file)

    if header != None:
        writer.writerow(header)

    for row in rows:
        writer.writerow(row)

    file.close()


def read_csv_file(file_name, skip_header_row = True):

    if skip_header_row:
        data = pandas.read_csv(file_name)
    else:
        data = pandas.read_csv(file_name, header=None)

    return data.values
