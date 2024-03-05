import csv

def create_position_dict(sequence):
    """Create a dictionary with position numbers for each letter in the sequence."""
    position_dict = {}
    for i, letter in enumerate(sequence, start=1):
        position_dict[i] = letter
    return position_dict

def align_isomiR_to_pre_miRNA(pre_miRNA_sequence, isomiR_sequence, isomiR_type):
    """Align isomiR sequence to pre-miRNA sequence based on isomiR type."""
    if isomiR_type.endswith('e3'):
        start_position = 0
    elif isomiR_type.endswith('e2'):
        start_position = 1
    elif isomiR_type == 'e1':
        start_position = 2
    elif isomiR_type == '3p':
        start_position = 3
    else:
        raise ValueError("Invalid isomiR type")

    aligned_sequence = ' ' * start_position + isomiR_sequence
    aligned_sequence += ' ' * (len(pre_miRNA_sequence) - len(aligned_sequence))

    return aligned_sequence

def match_letters(pre_miRNA_sequence, aligned_isomiR_sequence):
    """Match letters of aligned isomiR sequence with pre-miRNA sequence."""
    matched_letters = []
    for pre_letter, iso_letter in zip(pre_miRNA_sequence, aligned_isomiR_sequence):
        pre_letter = pre_letter.lower()
        iso_letter = iso_letter.lower()
        if iso_letter == ' ':
            matched_letters.append((' ', ' '))  # No match
        elif iso_letter == pre_letter:
            matched_letters.append((iso_letter, '+'))  # Match
        else:
            matched_letters.append((iso_letter, '-'))  # Mismatch
    return matched_letters

# Read input files
pre_miRNA_file = 'pre_miRNA.csv'
isomiR_file = 'isomiR_5_Ad.csv'

pre_miRNA_data = {}
with open(pre_miRNA_file, 'r') as f:
    reader = csv.reader(f)
    next(reader)  # Skip header
    for row in reader:
        name, sequence = row
        pre_miRNA_data[name] = sequence

isomiR_data = []
with open(isomiR_file, 'r') as f:
    reader = csv.reader(f)
    next(reader)  # Skip header
    for row in reader:
        isomiR_data.append(row)

# Process data and write output
output_file = 'output_5_Ad.csv'
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)

    # Write header
    writer.writerow(['Name', 'Pre-miRNA Sequence'] + [str(i) for i in range(1, len(max(pre_miRNA_data.values(), key=len)) + 1)])

    for name, pre_miRNA_sequence in pre_miRNA_data.items():
        writer.writerow([name, pre_miRNA_sequence] + list(pre_miRNA_sequence))

        for isomiR_name, isomiR_sequence, isomiR_type in isomiR_data:
            if isomiR_name == name:
                aligned_sequence = align_isomiR_to_pre_miRNA(pre_miRNA_sequence, isomiR_sequence, isomiR_type)
                matched_letters = match_letters(pre_miRNA_sequence, aligned_sequence)
                writer.writerow([isomiR_name, aligned_sequence] + list(matched_letters))
