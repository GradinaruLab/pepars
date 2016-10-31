import csv

def write_csv_file(file_name, header, rows):

    file = open(file_name, 'w')
    writer = csv.writer(file)
    writer.writerow(header)

    for row in rows:
        writer.writerow(row)

    file.close()

def read_csv_file(file_name, skip_header_row = True):

    file = open(file_name, 'r')
    reader = csv.reader(file, delimiter=',')

    data = []

    if skip_header_row:
        next(reader, None)

    count = 0

    for row in reader:
        column_iterator = iter(row)

        data_row = []
        for column in column_iterator:
            data_row.append(column)

        data.append(data_row)

        count += 1
        print(str(count))


    return data