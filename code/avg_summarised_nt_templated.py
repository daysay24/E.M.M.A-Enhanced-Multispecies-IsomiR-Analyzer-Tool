import pandas as pd 
import os
import sys
from colorama import Fore, Style, init
init(autoreset=True)

def get_avg(r, rep_cols):
    total_count = 0

    for rep_col in rep_cols: 
        total_count += r[rep_col]   

    return total_count / len(rep_cols)

def run(
    path_summarised_nt_alignment_output_folder,
    path_summarised_templated_alignment_output_folder,
    path_summarised_templated_alignment_all_output_folder,
    path_avg_summarised_nt_alignment_output_folder,
    path_avg_summarised_templated_alignment_output_folder,
    path_avg_summarised_templated_alignment_all_output_folder):
    print(Fore.MAGENTA + "\nAveraging statistics for different types of variation across multiple replicates ...")

    # Loop through each group
    for input_path, output_path in zip([path_summarised_nt_alignment_output_folder, path_summarised_templated_alignment_output_folder, path_summarised_templated_alignment_all_output_folder], [path_avg_summarised_nt_alignment_output_folder, path_avg_summarised_templated_alignment_output_folder, path_avg_summarised_templated_alignment_all_output_folder]):
        # Get group folders
        group_folders = os.listdir(input_path)
        # Create folder if not exists 
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        # Loop through each group  
        for group in group_folders:
            # Get the list of replicate files
            rep_files = os.listdir(f'{input_path}/{group}')
            # Create a dataframe that store replicates within the same group 
            group_df = pd.DataFrame()
            # Loop through each replicate file 
            for rep_file in rep_files:
                # Get replicate name 
                rep_name = rep_file.split('.')[0]
                # Read the replicate file 
                rep_df = pd.read_csv(f'{input_path}/{group}/{rep_file}', dtype={'Position': 'str'})
                # Key columns 
                key_cols = set(rep_df.columns) - {'value'}
                # Rename value column to replicate name
                rep_df = rep_df.rename(columns={'value': rep_name})
                # Check if the group_df is empty. If yes, group_df is set to be the first replicate 
                if group_df.empty:
                    group_df = rep_df
                # If not, merge that replicate to the current group_df 
                else:
                    group_df = pd.merge(group_df, rep_df, on=list(key_cols), how='outer')
                    group_df = group_df.fillna(0)

            # Get list of replicate columns 
            rep_cols = set(group_df.columns) - key_cols
            # Calculate the average value across all replicates
            group_df['count'] = group_df.apply(lambda r: get_avg(r, rep_cols), axis = 1)
            # Select subset of important columns  
            group_df = group_df[list(key_cols)+['count']]
            # Export to csv file 
            group_df.to_csv(f'{output_path}/{group}.csv', index=False)
                
