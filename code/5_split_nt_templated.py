import pandas as pd 
import os

def split_nt_templated(input_file, output_file, type):
    # Read input file
    templated_nt = pd.read_csv(input_file, low_memory=False)
    # Replace NA with ''
    templated_nt = templated_nt.fillna("(' ', ' ')")
    # Get position columns 
    cols = list(templated_nt.columns)
    del cols[0:3]
    # Loop over all rows 
    for index, r in templated_nt.iterrows():
        if r['is_pre'] == False:            
            for col in cols:
                vals = r[col]
                vals = vals.replace('(', '').replace(')', '').replace(' ', '').split(',')
                first_val = vals[0]
                second_val = vals[1]
                if first_val == '' and second_val == '':
                    templated_nt.loc[index, col] = ''
                elif type == 'templated':
                    templated_nt.loc[index, col] = second_val
                else: 
                    templated_nt.loc[index, col] = first_val  
    templated_nt.to_csv(output_file, index=False)

# Path to nt templated alignment files 
path_nt_templated_alignment_output_folder = '../data/4_nt_templated_alignment'
# Path to nt alignment output files 
path_nt_alignment_output_folder = '../data/5_nt_alignment'
# Path to templated alignment output files 
path_templated_alignment_output_folder = '../data/5_templated_alignment'
# List of group folders (e.g NEJ, JUV, AD)
group_folders = os.listdir(path_nt_templated_alignment_output_folder)
# Loop through each group
for group in group_folders:
    # Get the list of replicate files
    rep_files = os.listdir(f'{path_nt_templated_alignment_output_folder}/{group}')

    if not os.path.exists(f'{path_nt_alignment_output_folder}/{group}'):
        os.makedirs(f'{path_nt_alignment_output_folder}/{group}')

    if not os.path.exists(f'{path_templated_alignment_output_folder}/{group}'):
        os.makedirs(f'{path_templated_alignment_output_folder}/{group}')

    # Loop through each replicate file 
    for rep_file in rep_files:
        # Generate the nt alignment from nt templated file 
        split_nt_templated(f'{path_nt_templated_alignment_output_folder}/{group}/{rep_file}', f'{path_nt_alignment_output_folder}/{group}/{rep_file}', 'nt')
        # Generate the nt alignment from nt templated file 
        split_nt_templated(f'{path_nt_templated_alignment_output_folder}/{group}/{rep_file}', f'{path_templated_alignment_output_folder}/{group}/{rep_file}', 'templated')