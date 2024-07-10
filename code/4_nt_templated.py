import csv
import pandas as pd
import os

def align_isomiR_to_pre_miRNA(max_nt_diff_5p, nt_diff_5p, pre_seq, isomiR_seq):
    """Align isomiR sequence to pre-miRNA sequence based on isomiR type."""
    start_position = max_nt_diff_5p - nt_diff_5p
    aligned_seq = ' ' * start_position + isomiR_seq # start_position should be named gaps in pre-cursor
    aligned_seq += ' ' * (len(pre_seq) - len(aligned_seq))

    return aligned_seq

def match_letters(pri_seq, aligned_isomiR_seq):
    """Match letters of aligned isomiR sequence with pri-miRNA sequence."""
    matched_letters = []
    for pre_letter, iso_letter in zip(pri_seq, aligned_isomiR_seq):
        pre_letter = pre_letter.lower()
        iso_letter = iso_letter.lower()
        if iso_letter == ' ':
            matched_letters.append("(' ', ' ')")  # No match
        elif iso_letter == pre_letter:
            matched_letters.append(f"({iso_letter}, +)")  # Match
        else:
            matched_letters.append(f"({iso_letter}, -)")  # Mismatch
    return matched_letters

# Path to all summarised isomiR files 
path_summarised_output_folder = '../data/1_summarised_isomiRs'
# Path to extended precursor sequences file 
path_precursors_output_folder = '../data/3_precursors'
# Path to nt and templated alignment output files 
path_nt_templated_alignment_output_folder = '../data/4_nt_templated_alignment'
# List of group folders (e.g NEJ, JUV, AD)
group_folders = os.listdir(path_summarised_output_folder)
# Get precursor file 
precursor_output_file = [file for file in os.listdir(path_precursors_output_folder) if '.csv' in file][0]
# Read the extended precursor data  
extended_precursors = pd.read_csv(f'{path_precursors_output_folder}/{precursor_output_file}')
# Get max nt difference at 5p 
max_nt_diff_5p = int(precursor_output_file.split('_')[0])
# Loop through each group
for group in group_folders:
    # Get the list of replicate files
    rep_files = os.listdir(f'{path_summarised_output_folder}/{group}')

    if not os.path.exists(f'{path_nt_templated_alignment_output_folder}/{group}'):
        os.makedirs(f'{path_nt_templated_alignment_output_folder}/{group}')

    # Loop through each replicate file 
    for rep_file in rep_files:
        # Read the summarised isomiR-SEA output
        rep_df = pd.read_csv(f'{path_summarised_output_folder}/{group}/{rep_file}', encoding='latin-1')
        # Rename Pre-miRNA Name to mir_name
        rep_df = rep_df.rename(columns={'mirna_name': 'mir_name'})
        # Merge isomiR-SEA with extended_precursors to get the extended precursor sequence for each isomiR
        rep_df = rep_df.merge(extended_precursors, how='inner', on='mir_name')
        # Get rep name 
        rep_name = rep_file.split('.')[0]

        with open(f'{path_nt_templated_alignment_output_folder}/{group}/{rep_name}.csv', 'w+', newline='') as f:
            writer = csv.writer(f)
            # calculat max length of extended precursor
            max_extended_precursor_len = max([len(extended_precursor_seq) for extended_precursor_seq in list(extended_precursors['extended_precursor_seq'])])
            # Write header
            writer.writerow(['name', 'pre_seq', 'is_pre'] + [str(i) for i in range(1, max_extended_precursor_len + 1)])
            # Group isomiRs by mirna name and loop over each group 
            grouped_mir_name = rep_df.groupby('mir_name')
            for mir_name, mir_group in grouped_mir_name:
                # get the first record of mir_group 
                first_r = mir_group.iloc[0]
                pre_seq = first_r['extended_precursor_seq']
                # write the pri 
                writer.writerow([mir_name, pre_seq, True] + list(pre_seq))
                for _, r in mir_group.iterrows(): 
                    aligned_seq = align_isomiR_to_pre_miRNA(max_nt_diff_5p, r['5p_nt_diff'], pre_seq, r['tag_sequence'])
                    matched_letters = match_letters(pre_seq, aligned_seq)
                    writer.writerow([mir_name, aligned_seq, False] + list(matched_letters))


            
 