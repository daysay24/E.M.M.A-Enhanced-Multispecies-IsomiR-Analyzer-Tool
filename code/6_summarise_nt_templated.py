import pandas as pd 
from collections import Counter
import os

def get_extension_values(r, max_nt_diff_5p, max_nt_diff_3p):
    values = []
    pre_len = len(r['pre_seq'])
    for i in range(max_nt_diff_5p, 0, -1):
        values.append(r[str(i)])
    for i in range(pre_len - max_nt_diff_3p + 1, pre_len + 1):
        values.append(r[str(i)])
    return pd.Series(values)

def summarise_nt_alignment(path_nt_alignment_file, path_summarised_nt_alignment_file, max_nt_diff_5p, max_nt_diff_3p):
    # Read the nt alignment replicate file 
    nt_alignment = pd.read_csv(path_nt_alignment_file)

    # create a list of columns for extension positions at 5p
    extension_5p_cols = [ f"5'e{i + 1}" for i in range(max_nt_diff_5p)]
    # create a list of columns for extension positions at 3p
    extension_3p_cols = [ f"3'e{i + 1}" for i in range(max_nt_diff_3p)]
    # extension cols 
    extension_cols = extension_5p_cols + extension_3p_cols

    # Remove precursor rows 
    nt_alignment = nt_alignment[nt_alignment['is_pre'] == False]
    # For each row, get values for 5e1, 5e2, 5e3,..., 3e1, 3e2,... if that is an miRNA/isomiR 
    nt_alignment[extension_cols] = nt_alignment.apply(lambda r: get_extension_values(r, max_nt_diff_5p, max_nt_diff_3p), axis=1)
    # Summarise for nt 
    nt_summary = pd.DataFrame(columns=['position', 'nucleotide', 'value'])
    for col in extension_cols: 
        freq_counts = Counter(list(nt_alignment[col]))
        for nt in ['a', 'u', 'c', 'g']:
            nt_value = freq_counts[nt] if nt in freq_counts else 0
            nt_summary.loc[len(nt_summary.index)] = [col, nt, nt_value]
    nt_summary.to_csv(path_summarised_nt_alignment_file, index = False)

def summarise_templated_alignment(path_templated_alignment_file, path_summarised_templated_alignment_file, max_nt_diff_5p, max_nt_diff_3p):
    # Read the templated alignment replicate file 
    templated_alignment = pd.read_csv(path_templated_alignment_file)
    # create a list of columns for extension positions at 5p
    extension_5p_cols = [ f"5'e{i + 1}" for i in range(max_nt_diff_5p)]
    # create a list of columns for extension positions at 3p
    extension_3p_cols = [ f"3'e{i + 1}" for i in range(max_nt_diff_3p)]
    # extension cols 
    extension_cols = extension_5p_cols + extension_3p_cols

    # Remove precursor rows 
    templated_alignment = templated_alignment[templated_alignment['is_pre'] == False]
    # For each row, get values for 5e1, 5e2, 5e3,..., 3e1, 3e2,... if that is an miRNA/isomiR 
    templated_alignment[extension_cols] = templated_alignment.apply(lambda r: get_extension_values(r, max_nt_diff_5p, max_nt_diff_3p), axis=1)
    # Summarise for templated 
    templated_summary = pd.DataFrame(columns=['position', 'templated', 'value'])
    for col in extension_cols: 
        freq_counts = Counter(list(templated_alignment[col]))
        templated_value = freq_counts['+'] if '+' in freq_counts else 0
        untemplated_value = freq_counts['-'] if '+' in freq_counts else 0
        templated_summary.loc[len(templated_summary.index)] = [col, 'Templated', templated_value]
        templated_summary.loc[len(templated_summary.index)] = [col, 'Untemplated', untemplated_value]
        templated_summary.to_csv(path_summarised_templated_alignment_file, index = False)

def summarise_templated_alignment_all(path_templated_alignment_file, path_summarised_templated_alignment_all_file, max_nt_diff_5p):
    # Read the templated alignment replicate file 
    templated_alignment = pd.read_csv(path_templated_alignment_file)
    # Remove precursor rows 
    templated_alignment = templated_alignment[templated_alignment['is_pre'] == False]
    # Retain extended or truncated isomiRs 
    templated_alignment = templated_alignment[templated_alignment['extended_or_truncated'].isin(['extended', 'truncated'])]
    # Summarise for templated 
    templated_summary = pd.DataFrame(columns=['type', 'position', 'templated', 'value'])
    # Group by extended / truncated 
    grouped_extended_or_truncated = templated_alignment.groupby('extended_or_truncated')
    # Position columns 
    cols = list(templated_alignment.columns)
    del cols[0:4]
    for extended_or_truncated, extended_or_truncated_group in grouped_extended_or_truncated:
        # Loop through each position 
        for col in cols: 
            freq_counts = Counter(list(extended_or_truncated_group[col]))
            templated_value = freq_counts['+'] if '+' in freq_counts else 0
            untemplated_value = freq_counts['-'] if '+' in freq_counts else 0
            col = int(col)
            if extended_or_truncated == 'extended':
                col = f"5'e{max_nt_diff_5p - col + 1}" if col <= max_nt_diff_5p else col - max_nt_diff_5p
            else: 
                if col <= max_nt_diff_5p:
                    continue
                else: 
                    col = col - max_nt_diff_5p
                
            templated_summary.loc[len(templated_summary.index)] = [extended_or_truncated, col, 'Templated', templated_value]
            templated_summary.loc[len(templated_summary.index)] = [extended_or_truncated, col, 'Untemplated', untemplated_value]
    templated_summary.to_csv(path_summarised_templated_alignment_all_file, index = False)

# Path to nt alignment output files 
path_nt_alignment_output_folder = '../data/5_nt_alignment'
# Path to templated alignment output files 
path_templated_alignment_output_folder = '../data/5_templated_alignment'
# Path to summarised nt alignment output files 
path_summarised_nt_alignment_output_folder = '../data/6_summarised_nt_alignment'
# Path to summarised templated alignment output files 
path_summarised_templated_alignment_output_folder = '../data/6_summarised_templated_alignment'
# Path to summarised templated alignment all positions output files 
path_summarised_templated_alignment_all_output_folder = '../data/6_summarised_templated_alignment_all'
# List of group folders (e.g NEJ, JUV, AD)
nt_group_folders = os.listdir(path_nt_alignment_output_folder)
templated_group_folders = os.listdir(path_templated_alignment_output_folder)
# Path to extended precursor sequences file 
path_precursors_output_folder = '../data/3_precursors'
# Get precursor file 
precursor_output_file = [file for file in os.listdir(path_precursors_output_folder) if '.csv' in file][0]
# Get max nt difference at 5p 
max_nt_diff_5p, max_nt_diff_3p = int(precursor_output_file.split('_')[0]), int(precursor_output_file.split('_')[1])

for nt_group, templated_group in zip(nt_group_folders, templated_group_folders):
    # Get the list of nt alignment replicate files
    rep_nt_files = os.listdir(f'{path_nt_alignment_output_folder}/{nt_group}')
    # Get the list of templated alignment replicate files
    rep_templated_files = os.listdir(f'{path_templated_alignment_output_folder}/{templated_group}')

    if not os.path.exists(f'{path_summarised_nt_alignment_output_folder}/{nt_group}'):
        os.makedirs(f'{path_summarised_nt_alignment_output_folder}/{nt_group}')

    if not os.path.exists(f'{path_summarised_templated_alignment_output_folder}/{templated_group}'):
        os.makedirs(f'{path_summarised_templated_alignment_output_folder}/{templated_group}')

    if not os.path.exists(f'{path_summarised_templated_alignment_all_output_folder}/{templated_group}'):
        os.makedirs(f'{path_summarised_templated_alignment_all_output_folder}/{templated_group}')

    # Loop through each replicate file 
    for nt_rep_file, templated_rep_file in zip(rep_nt_files, rep_templated_files):
        summarise_nt_alignment(f'{path_nt_alignment_output_folder}/{nt_group}/{nt_rep_file}', 
                               f'{path_summarised_nt_alignment_output_folder}/{nt_group}/{nt_rep_file}',
                               max_nt_diff_5p,
                               max_nt_diff_3p)
        summarise_templated_alignment(f'{path_templated_alignment_output_folder}/{templated_group}/{templated_rep_file}', 
                               f'{path_summarised_templated_alignment_output_folder}/{templated_group}/{templated_rep_file}',
                               max_nt_diff_5p,
                               max_nt_diff_3p)
        summarise_templated_alignment_all(f'{path_templated_alignment_output_folder}/{templated_group}/{templated_rep_file}', 
                                f'{path_summarised_templated_alignment_all_output_folder}/{templated_group}/{templated_rep_file}',
                                max_nt_diff_5p)
                
                                                                                                