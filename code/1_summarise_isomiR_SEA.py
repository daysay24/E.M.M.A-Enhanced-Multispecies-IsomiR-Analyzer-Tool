# This code will summarise raw outputs from isomiR-SEA for all replicates

import pandas as pd 
import os


def snp_count(mir_seq, tag_seq, begin_ungapped_mirna, begin_ungapped_tag):
    snp_count = 0
    mir_seq_nts = list(mir_seq)
    tag_seq_nts = list(tag_seq)
    mir_seq_len = len(mir_seq_nts)
    tag_seq_len = len(tag_seq_nts)

    while begin_ungapped_mirna < mir_seq_len - 1 and begin_ungapped_tag < tag_seq_len - 1:
        if mir_seq_nts[begin_ungapped_mirna] != tag_seq_nts[begin_ungapped_tag]:
            snp_count += 1
        begin_ungapped_mirna += 1
        begin_ungapped_tag += 1
    
    return snp_count

def get_type_annotation(mir_name, nt_diff_5p, nt_snp, nt_diff_3p):    
    type = []
    annotation = []

    if nt_diff_5p == 0 and nt_snp == 0 and nt_diff_3p == 0:
        return 'mirna_exact', mir_name

    if nt_diff_5p > 0: 
        annotation.append(f"5'+{nt_diff_5p}")
        type.append('iso_5p')
    elif nt_diff_5p < 0: 
        annotation.append(f"5'{nt_diff_5p}")
        type.append('iso_5p')

    if nt_snp > 1: 
        annotation.append(f"snp+{nt_snp}")
        type.append('iso_multi_snp')  
    elif nt_snp == 1: 
        annotation.append(f"snp+{nt_snp}")
        type.append('iso_snp')  
        
    if nt_diff_3p > 0: 
        annotation.append(f"3'+{nt_diff_3p}")
        type.append('iso_3p')
    elif nt_diff_3p < 0: 
        annotation.append(f"3'{nt_diff_3p}")
        type.append('iso_3p')

    type_str = '-'.join(type) if len(type) > 1 else '-'.join(type) + "_only"

    return type_str, mir_name + '(' + '|'.join(annotation) + ')'

def add_columns(row):
    # 5p nt diff 
    nt_diff_5p = row['begin_ungapped_tag'] - row['begin_ungapped_mirna']
    # 3p nt diff 
    nt_diff_3p = - row['mir_tag_size_diff'] - nt_diff_5p
    # snps
    nt_snp = snp_count(row['mirna_seq'], row['tag_sequence'], row['begin_ungapped_mirna'], row['begin_ungapped_tag'])
    # type and annotation
    type, annotation = get_type_annotation(row['mirna_name'], nt_diff_5p, nt_snp, nt_diff_3p)

    # add columnns
    return pd.Series([nt_diff_5p, nt_snp, nt_diff_3p, type, annotation])

# Path to the raw outputs folder 
path_raw_output_folder = '../data/0_isomiR-SEA_isomiRs'
path_summarised_output_folder = '../data/1_summarised_isomiRs'
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
        # Get the mir name  
        isomiR_SEA_output['mirna_name'] = isomiR_SEA_output['mirna_name'].apply(lambda r: r.replace('>', '').split(' ')[0])
        # Add 5p_nt_diff, snp_nt, 3p_nt_diff, type, annotation
        isomiR_SEA_output[['5p_nt_diff', 'snp_nt', '3p_nt_diff', 'type', 'annotation']] = isomiR_SEA_output.apply(lambda r: add_columns(r), axis = 1)
        # Create folder if not exist 
        if not os.path.exists(f'{path_summarised_output_folder}/{group}'):
            os.makedirs(f'{path_summarised_output_folder}/{group}')
        isomiR_SEA_output.to_csv(f'{path_summarised_output_folder}/{group}/{rep_file}', index=False)
