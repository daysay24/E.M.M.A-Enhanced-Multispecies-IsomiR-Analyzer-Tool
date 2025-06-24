import csv
import pandas as pd
import os
import sys
from colorama import Fore, Style, init
init(autoreset=True)

def align_isomiR_to_pre_miRNA(max_nt_diff_5p, nt_diff_5p, pre_seq, isomiR_seq):
    """Align an isomiR sequence to an extended precursor sequence. 
    The function returns the new isomiR sequence having the same length with the extended precursor sequence, in which unaligned positions are filled with spaces.

    Parameters
    ----------
    max_nt_diff_5p : int
        The maximum number of nucleotide difference at 5' end across all isomiRs.
    nt_diff_5p : int 
        The number of nucleotide added to / trimmed from the 5' of the canonical. 
    pre_seq : str 
        The extended precursor sequence of a miRNA, which is extracted by adding max_nt_diff_5p and max_nt_diff_3p nucleotids to 5' and 3' ends of that miRNA. 
    isomiR_seq : str 
        The isomiR sequence. 

    Example
    -----------
    ```
    Extended precursor sequence : "AUGCUAUCCCGCUAAUGCUAUCCCGCU"
    
    isomiR :                          "UAUCCCGCUAAUGCUAU"      

    Output:                       "    UAUCCCGCUAAUGCUAU     "
    ```
    
    Returns
    ------
    str 
        An isomiR sequence aligned to the extended precursor sequence. 
    """
    start_position = max_nt_diff_5p - nt_diff_5p
    aligned_seq = ' ' * start_position + isomiR_seq # start_position should be named gaps in pre-cursor
    aligned_seq += ' ' * (len(pre_seq) - len(aligned_seq))

    return aligned_seq

def match_letters(pre_seq, aligned_isomiR_seq):
    """Convert each nucleotide in isomiR to (<nucleotide>, <matching symbol>) format by comparing the isomiR with the extended precursor. 
    
    At each position of an isomiR: 
    - The nucleotide of the isomiR is the same as that of the extended precursor sequence, matching symbol is +.
    - The nucleotide of the isomiR is different from that of the extended precursor sequence, matching symbol is -.
    
    Parameters
    ----------
    pre_seq : str
        The extended precursor sequence of a miRNA, which is extracted by adding max_nt_diff_5p and max_nt_diff_3p nucleotides to 5' and 3' ends of that miRNA. 
    aligned_isomiR_seq : str
        The output of align_isomiR_to_pre_miRNA().

    Example 
    -----------
    ```
    Extended precursor sequence : "AUGCUAUCCUGCUGUCCCG"
    
    isomiR :                      "    GAUCCUGCUAU    "   

    Output: [(' ', ' '),  (' ', ' '), (' ', ' '), (' ', ' '), ('g', '-'), ('a', '+'), ('u', '+'), ('c', '+'), ('c', '+'), ('u', '+'), ('g', '+'), ('c', '+'), ('u', '+'), ('g', '-'), ('u', '+'), (' ', ' '),  (' ', ' '), (' ', ' '), (' ', ' ')]  
    ``` 

    Returns
    -------
    List 
        A list of nucleotide details in (<nucleotide>, <matching symbol>) format at each position of the isomiR. 
    """
    matched_letters = []
    for pre_letter, iso_letter in zip(pre_seq, aligned_isomiR_seq):
        pre_letter = pre_letter.lower()
        iso_letter = iso_letter.lower()
        if iso_letter == ' ':
            matched_letters.append("(' ', ' ')")  # No match
        elif iso_letter == pre_letter:
            matched_letters.append(f"({iso_letter}, +)")  # Match
        else:
            matched_letters.append(f"({iso_letter}, -)")  # Mismatch
    return matched_letters

def extended_or_truncated(nt_5p_diff, nt_3p_diff):
    """
    Classify an isomiR into truncated, extended, or neither. 

    Parameters
    ----------
    nt_5p_diff : int
        The number of nucleotide added to / trimmed from the 5' of the canonical. 
    nt_3p_diff : int
        The number of nucleotide added to / trimmed from the 3' of the canonical. 

    - If at least one of the end is truncated, the output is 'truncated'.
    - If at least one of the end is extended, the output is 'extended'.
    - Otherwise, the output is ''.
    Returns
    -------
    str: 
        Type of an isomiR based on alterations at 5' and 3'. 
    
    """
    if (nt_5p_diff < 0 and nt_3p_diff == 0) or (nt_5p_diff == 0 and nt_3p_diff < 0) or (nt_5p_diff < 0 and nt_3p_diff < 0):
    # if (nt_5p_diff < 0 or nt_3p_diff < 0):
        return 'truncated'
    elif (nt_5p_diff > 0 and nt_3p_diff == 0) or (nt_5p_diff == 0 and nt_3p_diff > 0) or (nt_5p_diff > 0 and nt_3p_diff > 0):
    # elif (nt_5p_diff > 0 or nt_3p_diff > 0):
        return 'extended'
    else: 
        return ''

def run(path_summarised_output_folder, path_precursors_output_folder, path_nt_templated_alignment_output_folder):
    print(Fore.MAGENTA + "\nComparing nucleotide at each position of isomiRs ...")

    # List of group folders 
    group_folders = os.listdir(path_summarised_output_folder)
    # Get precursor file 
    precursor_output_file = [file for file in os.listdir(path_precursors_output_folder) if '.csv' in file][0]
    # Read the extended precursor data  
    extended_precursors = pd.read_csv(f'{path_precursors_output_folder}/{precursor_output_file}')
    # Get max nt difference at 5p 
    max_nt_diff_5p = int(precursor_output_file.split('_')[0])
    # Loop through each group folder
    for group in group_folders:
        # Get the list of replicate files
        rep_files = os.listdir(f'{path_summarised_output_folder}/{group}')

        if not os.path.exists(f'{path_nt_templated_alignment_output_folder}/{group}'):
            os.makedirs(f'{path_nt_templated_alignment_output_folder}/{group}')

        # Loop through each replicate file 
        for rep_file in rep_files:
            # Read the replicate file
            rep_df = pd.read_csv(f'{path_summarised_output_folder}/{group}/{rep_file}', encoding='latin-1')
            # Rename mirna_name to mir_name
            rep_df = rep_df.rename(columns={'mirna_name': 'mir_name'})
            # Merge with extended_precursors to get the extended precursor sequence for each isomiR
            rep_df = rep_df.merge(extended_precursors, how='inner', on='mir_name')
            # Get replicate name 
            rep_name = rep_file.split('.')[0]

            with open(f'{path_nt_templated_alignment_output_folder}/{group}/{rep_name}.csv', 'w+', newline='') as f:
                writer = csv.writer(f)
                # Calculate max length of extended precursor
                max_extended_precursor_len = max([len(extended_precursor_seq) for extended_precursor_seq in list(extended_precursors['extended_precursor_seq'])])
                # Write file header
                writer.writerow(['name', 'pre_seq', 'is_pre', 'extended_or_truncated'] + [str(i) for i in range(1, max_extended_precursor_len + 1)])
                # Group isomiRs by mirna name and loop over each group 
                grouped_mir_name = rep_df.groupby('mir_name')
                for mir_name, mir_group in grouped_mir_name:
                    # Get the first record of mir_group 
                    first_r = mir_group.iloc[0]
                    pre_seq = first_r['extended_precursor_seq']
                    # Write a row for the extended precursor for the miRNA.
                    writer.writerow([mir_name, pre_seq, True, ''] + list(pre_seq))
                    for _, r in mir_group.iterrows(): 
                        aligned_seq = align_isomiR_to_pre_miRNA(max_nt_diff_5p, r['5p_nt_diff'], pre_seq, r['tag_sequence'])
                        matched_letters = match_letters(pre_seq, aligned_seq)
                        writer.writerow([mir_name, aligned_seq, False, extended_or_truncated(r['5p_nt_diff'], r['3p_nt_diff'])] + list(matched_letters))


                
    