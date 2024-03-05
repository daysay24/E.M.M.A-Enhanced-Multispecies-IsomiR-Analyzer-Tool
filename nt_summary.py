import csv

def process_txt(input_file, output_file):
    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        data = list(reader)

    # Process the data to remove "(' ', '-')" and format accordingly
    for row_index, row in enumerate(data[2:], start=2):
        for i in range(2, len(row)):
            if row[i] == "(' ', '-')":
                row[i] = ' '
            else:
                cleaned_data = row[i].replace('(', '').replace(')', '').replace("'", "").replace(",", "").replace("+", "").split()
                if cleaned_data:
                    row[i] = cleaned_data[0]
                else:
                    print(f"Warning: Insufficient data in row {row_index}, element {i}")

    # Write processed data to the output file
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerows(data)

input_file = "summary_input_5.txt"
output_file = "output_nt_summary_5.txt"
process_txt(input_file, output_file)
