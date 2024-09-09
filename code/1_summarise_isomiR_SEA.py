import pandas as pd 
import os
import sys

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
    
    while mirna_aligned_start < mir_seq_len and tag_aligned_start < tag_seq_len:
        if mir_seq_nts[mirna_aligned_start] != tag_seq_nts[tag_aligned_start]:
            snp_count += 1
        mirna_aligned_start += 1
        tag_aligned_start += 1
    
    return snp_count

def get_type_name(mir_name: str, nt_diff_5p: int, nt_snp: int, nt_diff_3p: int):    
    """Get the variant type and name of an isomiR.

    Parameters
    ----------
    mir_name: str
        Name of its canonical miRNA. 
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
        name.append(f"5'+{nt_diff_5p}")
        type.append('iso_5p')
    elif nt_diff_5p < 0: 
        name.append(f"5'{nt_diff_5p}")
        type.append('iso_5p')

    if nt_snp > 1: 
        name.append(f"snp+{nt_snp}")
        type.append('iso_multi_snp')  
    elif nt_snp == 1: 
        name.append(f"snp+{nt_snp}")
        type.append('iso_snp')  
        
    if nt_diff_3p > 0: 
        name.append(f"3'+{nt_diff_3p}")
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
    nt_snp = snp_count(row['mirna_seq'], row['tag_sequence'], row['begin_ungapped_mirna'], row['begin_ungapped_tag'])
    type, name = get_type_name(row['mirna_name'], nt_diff_5p, nt_snp, nt_diff_3p)

    return pd.Series([nt_diff_5p, nt_snp, nt_diff_3p, type, name])

# Path to the raw isomiR outputs folder 
path_raw_output_folder = sys.argv[1]
# Path to the summarised isomiR outputs folder
path_summarised_output_folder = sys.argv[2]

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
        # Add 5p_nt_diff, snp_nt, 3p_nt_diff, type, annotation columns 
        isomiR_SEA_output[['5p_nt_diff', 'snp_nt', '3p_nt_diff', 'type', 'annotation']] = isomiR_SEA_output.apply(lambda r: add_columns(r), axis = 1)
        
        # Create folder if not exist 
        if not os.path.exists(f'{path_summarised_output_folder}/{group}'):
            os.makedirs(f'{path_summarised_output_folder}/{group}')
        isomiR_SEA_output.to_csv(f'{path_summarised_output_folder}/{group}/{rep_file}', index=False)