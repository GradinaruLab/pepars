import csv

def write_csv_file(file_name, header, rows):

    file = open(file_name, 'w')
    writer = csv.writer(file)
    writer.writerow(header)

    for row in rows:
        writer.writerow(row)

    file.close()