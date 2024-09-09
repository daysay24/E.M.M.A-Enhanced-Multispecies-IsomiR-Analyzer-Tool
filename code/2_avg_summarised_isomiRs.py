import pandas as pd 
import os 
import sys

def get_grouped_type(type: str):
    """Categorise 12 variant types into 5 big groups 
    
    Parameters
    ----------
    type: str 
        Variant type of an isomiR. 
        - 'iso_3p_only' -> 3'isomiR
        - 'iso_5p_only' -> 5'isomiR
        - 'iso_5p-iso_snp-iso_3p', 'iso_5p-iso_multi_snp-iso_3p', 'iso_5p-iso_3p' -> 'Both end isomiR'
        - 'mirna_exact' -> Canonical
        - The rest -> Others
    Returns
    -------
    str 
        The group of variant type of an isomiR. 

    """
    if type == 'iso_3p_only':
        return "3'isomiR"
    elif type == 'iso_5p_only':
        return "5'isomiR"
    elif type in ['iso_5p-iso_snp-iso_3p', 'iso_5p-iso_multi_snp-iso_3p', 'iso_5p-iso_3p']:
        return 'Both end isomiR'
    elif type == 'mirna_exact':
        return 'Canonical'
    else: 
        return 'Others'

def get_type_nt(type: str, nt_diff_5p: int, nt_diff_3p: int):
    """Get the combination of variant type and the number of nucleotide differences of an iso_5p_only or iso_3p_only isomiR. 
    
    Paramaters
    ----------
    type : str 
        The variant type of the isomiR. 
    nt_diff_5p: int
        The number of nucleotide added to / trimmed from the 5' of the canonical. 
    nt_diff_3p: int
        The number of nucleotide added to / trimmed from the 3' of the canonical. 

    Returns
    -------
    str
        The combination of variant type and the number of nucleotide differences of an isomiR.

    """
    if type in ['iso_3p_only', 'iso_5p_only']:
        side, nt = '', ''
        if nt_diff_5p != 0:  
            side, nt = 5, nt_diff_5p
        else: 
            side, nt = 3, nt_diff_3p 
        type_nt = f"{side}'t{nt}" if nt < 0 else f"{side}'e{nt}"
        return type_nt
    else: 
        return ''
    
def get_avg(r: pd.Series, rep_cols: set):
    """Calculate the average rpm and unique tag for each isomiR across multiple replicates

    Parameters
    ----------
    r : pandas.Series
        A 1D array that stores information of an isomiR. 
    rep_cols : set
        A set of columns that store rpm / unique tag count of an isomiR.  

    Returns
    -------
    pd.Series
        The average rpm and average unique tag count of an isomiR across multiple replicates. 
    """
    total_rpm = 0
    total_unique_tag = 0
    n_reps = len(rep_cols) / 2

    for rep_col in rep_cols: 
        if 'rpm' in rep_col:
            total_rpm += r[rep_col]
        else:
            total_unique_tag += r[rep_col]       

    return pd.Series([total_rpm / n_reps, total_unique_tag / n_reps])     

# Path to all summarised isomiR files 
path_summarised_output_folder = sys.argv[1]
# Path to the average rpm/unique tag isomiR files 
path_avg_replicate_output_folder = sys.argv[2]

# List of group folders 
group_folders = os.listdir(path_summarised_output_folder)
# Loop through each group
for group in group_folders:
    # Get the list of replicate files
    rep_files = os.listdir(f'{path_summarised_output_folder}/{group}')
    # Create a dataframe that store summarised isomiRs of all replicates within the same group 
    group_df = pd.DataFrame()
    # Create folder if not exists 
    if not os.path.exists(path_avg_replicate_output_folder):
        os.makedirs(path_avg_replicate_output_folder)

    # Loop through each replicate file 
    for rep_file in rep_files:
        # Get replicate name 
        rep_name = rep_file.split('.')[0]
        # Read the replicate file 
        rep_df = pd.read_csv(f'{path_summarised_output_folder}/{group}/{rep_file}')
        # Select a subset of important columns 
        rep_df = rep_df[['mirna_name', 'tag_sequence', 'type', '#count_tags', '5p_nt_diff', '3p_nt_diff']]
        # Normalise raw count and store in a new column named <replicate_name>_rpm
        sum_raw_count = sum(list(rep_df['#count_tags']))
        rep_df[f'{rep_name}_rpm'] = rep_df['#count_tags'] * 1000000 / sum_raw_count
        # Remove the #count_tags column 
        rep_df = rep_df.drop(columns=['#count_tags'])
        # Add unique tag column and set value to 1
        rep_df[f'{rep_name}_unique_tag'] = 1
        # Group variant types into 3p, 5p, both, canonical and others and save to a new column grouped_type
        rep_df['grouped_type'] = rep_df['type'].apply(lambda t: get_grouped_type(t))
        # Combine type and the number of nt differences and save to a new column type_nt 
        rep_df['type_nt'] = rep_df.apply(lambda r: get_type_nt(r['type'], r['5p_nt_diff'], r['3p_nt_diff']), axis = 1)
        # Check if the group_df is empty. If yes, group_df is set to be the summarised isomiRs of the first replicate 
        if group_df.empty:
            group_df = rep_df
        # if not, merge the summarised isomiRs of that replicated to the current group_df 
        else:
            group_df = pd.merge(group_df, rep_df, on=['mirna_name', 'tag_sequence', 'type', '5p_nt_diff', '3p_nt_diff', 'grouped_type', 'type_nt'], how='outer')
            group_df = group_df.fillna(0)

    # Get list of replicate columns that store rpm / unique tag count of isomiRs 
    rep_cols = set(group_df.columns) - {'mirna_name', 'tag_sequence', 'type', '5p_nt_diff', '3p_nt_diff', 'grouped_type', 'type_nt'}
    # Calculate the average rpm and unique tag for each isomiR and save to a new columns rpm, unique_tag
    group_df[['rpm', 'unique_tag']] = group_df.apply(lambda r: get_avg(r, rep_cols), axis = 1)
    # Select subset of important columns  
    group_df = group_df[['mirna_name', 'tag_sequence', 'grouped_type', 'type_nt', 'rpm', 'unique_tag']]
    # Export to csv file 
    group_df.to_csv(f'{path_avg_replicate_output_folder}/{group}.csv', index=False)