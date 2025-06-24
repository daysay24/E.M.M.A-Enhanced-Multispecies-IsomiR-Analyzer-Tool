import pandas as pd 
import os
import sys
from colorama import Fore, Style, init
init(autoreset=True)

def snp_count(mir_seq: str, tag_seq: str, begin_ungapped_mirna: int, begin_ungapped_tag: int):
    """Count the number of snp found in an isomiR. 
    
    Parameters 
    ----------
    mir_seq : str
        Sequence of its canonical miRNA.
    tag_seq: str 
        Sequence of the isomiR. 
    begin_ungapped_mirna: int
        The miRNA nucleotide position where the tag alignment begin.
    begin_ungapped_tag: int 
        The tag nucleotide position where the alignment with the selected miRNA begin.

    Returns
    -------
    int 
        The number of snps found in an isomiR. 
        The string that describes position and nt of snps. 
    """
    snp_count = 0
    mir_seq_nts = list(mir_seq)
    tag_seq_nts = list(tag_seq)
    mir_seq_len = len(mir_seq_nts)
    tag_seq_len = len(tag_seq_nts)

    # mirna position where mirna is first aligned to tag even if mismatch
    mirna_aligned_start = begin_ungapped_mirna - min(begin_ungapped_mirna, begin_ungapped_tag)
    # tag position where tag is first aligned to mirna even if mismatch
    tag_aligned_start = begin_ungapped_tag - min(begin_ungapped_mirna, begin_ungapped_tag)
    # snp description 
    snp_desc = ''
    
    while mirna_aligned_start < mir_seq_len and tag_aligned_start < tag_seq_len:
        if mir_seq_nts[mirna_aligned_start] != tag_seq_nts[tag_aligned_start]:
            snp_count += 1
            snp_desc += f"{mirna_aligned_start + 1}{tag_seq_nts[tag_aligned_start]}"

        mirna_aligned_start += 1
        tag_aligned_start += 1
    
    return snp_count, snp_desc

def get_type_name(mir_name: str, tag_seq: str, nt_diff_5p: int, nt_snp: int, snp_desc: str,  nt_diff_3p: int):    
    """Get the variant type and name of an isomiR.

    Parameters
    ----------
    mir_name: str
        Name of its canonical miRNA. 
    tag_seq: str 
        Sequence of the isomiR. 
    nt_diff_5p: int
        The number of nucleotide added to / trimmed from the 5' of the canonical. 
        - nt_diff_5p = 0 (no addition/trimming)
        - nt_diff_5p > 0 (addition)
        - nt_diff_5p < 0 (trimming)
    nt_snp: int
        The number of snp found in the isomiR.  
        - nt_snp = 0 (no snp)
        - nt_snp <> 0 (snp)
    nt_diff_3p: int
        The number of nucleotide added to / trimmed from the 3' of the canonical. 
        - nt_diff_3p = 0 (no addition/trimming)
        - nt_diff_3p > 0 (addition)
        - nt_diff_3p < 0 (trimming)

    Returns
    -------
    str, str
        Variant type and name of an isomiR.

    """
    type = []
    name = []

    if nt_diff_5p == 0 and nt_snp == 0 and nt_diff_3p == 0:
        return 'mirna_exact', mir_name

    if nt_diff_5p > 0: 
        name.append(f"5'+{nt_diff_5p}:{tag_seq[:nt_diff_5p]}")
        type.append('iso_5p')
    elif nt_diff_5p < 0: 
        name.append(f"5'{nt_diff_5p}")
        type.append('iso_5p')

    if nt_snp > 1: 
        name.append(f"snp+{nt_snp}:{snp_desc}")
        type.append('iso_multi_snp')  
    elif nt_snp == 1: 
        name.append(f"snp+{nt_snp}:{snp_desc}")
        type.append('iso_snp')  
        
    if nt_diff_3p > 0: 
        name.append(f"3'+{nt_diff_3p}:{tag_seq[-nt_diff_3p:]}")
        type.append('iso_3p')
    elif nt_diff_3p < 0: 
        name.append(f"3'{nt_diff_3p}")
        type.append('iso_3p')

    type_str = '-'.join(type) if len(type) > 1 else type[0] + "_only"

    return type_str, mir_name + '(' + '|'.join(name) + ')'

def add_columns(row: pd.Series):
    """Add extra information to isomiR-SEA output for downstream analysis. 

    Parameters
    ----------
    row: pandas.Series
        A 1D array that stores values of a row in the isomiR-SEA output table. Indexes correspond to column names, values correspond to cell values.  
    
    Returns
    -------
    pandas.Series
        A 1D array that stores extra information of a row in the isomiR-SEA output table. 5 extra values: 
        - nt_diff_5p: The number of nucleotide added to / trimmed from the 5' of the canonical. 
        - nt_snp: The number of snp found in the isomiR.  
        - nt_diff_3p: The number of nucleotide added to / trimmed from the 3' of the canonical. 
        - type: The variant type of the isomiR. 
        - name: The name of the isomiR. 
    """
    nt_diff_5p = row['begin_ungapped_tag'] - row['begin_ungapped_mirna']
    nt_diff_3p = - row['mir_tag_size_diff'] - nt_diff_5p
    nt_snp, snp_desc = snp_count(row['mirna_seq'], row['tag_sequence'], row['begin_ungapped_mirna'], row['begin_ungapped_tag'])
    type, name = get_type_name(row['mirna_name'], row['tag_sequence'],  nt_diff_5p, nt_snp, snp_desc, nt_diff_3p)

    return pd.Series([nt_diff_5p, nt_snp, nt_diff_3p, type, name])

def run(path_raw_output_folder, path_summarised_output_folder, read_count_threshold):
    print(Fore.MAGENTA + "\nUpdating outputs of isomiR-SEA by calculating 5', 3' and snp modification, naming isomiRs, categorizing isomiRs, ...")

    # List of all tag sequences and their sum of read counts across all samples e.g {'AACCCUGUAGACCCGAGUUUGG': 34, 'UGAAAGACGAUGGUAGUGAGAUG': 10, 'ACCCUUGUUCGACUGUGA': 8, ...}
    sum_read_counts = {}

    # List of all sample groups 
    group_folders = os.listdir(path_raw_output_folder)
    # Loop through each group
    for group in group_folders:
        # Get the list of replicate files
        rep_files = os.listdir(f'{path_raw_output_folder}/{group}')
        # Loop through each replicate of that group 
        for rep_file in rep_files:
            # Read isomiR-SEA raw output file of that replicate
            isomiR_SEA_output = pd.read_csv(f'{path_raw_output_folder}/{group}/{rep_file}', sep='\t', encoding='latin-1')
            # Change the datatype of the column 
            isomiR_SEA_output = isomiR_SEA_output.astype({
                'begin_ungapped_tag': int,
                'begin_ungapped_mirna': int,
                'mir_tag_size_diff': int
            })
            # Retain the mirna_name   
            isomiR_SEA_output['mirna_name'] = isomiR_SEA_output['mirna_name'].apply(lambda r: r.replace('>', '').split(' ')[0])
            # Exclude reads having N inside it
            isomiR_SEA_output = isomiR_SEA_output[isomiR_SEA_output['tag_sequence'].str.contains('N') == False]
            # Add 5p_nt_diff, snp_nt, 3p_nt_diff, type, annotation columns 
            isomiR_SEA_output[['5p_nt_diff', 'snp_nt', '3p_nt_diff', 'type', 'annotation']] = isomiR_SEA_output.apply(lambda r: add_columns(r), axis = 1)
            
            # Get list of tag sequences found so far 
            found_tag_sequences = list(sum_read_counts.keys())
            # Loop through each tag sequence
            for _, r in isomiR_SEA_output.iterrows(): 
                if r['tag_sequence'] in found_tag_sequences: 
                    sum_read_counts[r['tag_sequence']] = sum_read_counts[r['tag_sequence']] + r['#count_tags']
                else: 
                    sum_read_counts[r['tag_sequence']] = r['#count_tags']

            # Create folder if not exist 
            if not os.path.exists(f'{path_summarised_output_folder}/{group}'):
                os.makedirs(f'{path_summarised_output_folder}/{group}')
            isomiR_SEA_output.to_csv(f'{path_summarised_output_folder}/{group}/{rep_file}', index=False)

    # Get tag sequences having read counts >= read_count_threshold
    kept_tag_sequences = [k for k,v in sum_read_counts.items() if v >= read_count_threshold]
    # Loop through each group
    for group in group_folders:
        # Get the list of replicate files
        rep_files = os.listdir(f'{path_summarised_output_folder}/{group}')
        # Loop through each replicate of that group 
        for rep_file in rep_files:
            # Read summarised isomiR-SEA output file of that replicate
            isomiR_SEA_output = pd.read_csv(f'{path_summarised_output_folder}/{group}/{rep_file}')
            # Keep tag sequences having total read counts >= read_count_threshold 
            isomiR_SEA_output = isomiR_SEA_output[isomiR_SEA_output['tag_sequence'].isin(kept_tag_sequences)]
            isomiR_SEA_output.to_csv(f'{path_summarised_output_folder}/{group}/{rep_file}', index=False)