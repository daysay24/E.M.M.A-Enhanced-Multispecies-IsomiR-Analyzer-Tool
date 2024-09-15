import pandas as pd 
import os
import sys

def split_nt_templated(input_file, output_file, type):
    """Extract nucleotide or matching symbol from the nucleotide details (<nucleotide>, <matching symbol>) and store each in a seperate file. 

    Parameters 
    ----------
    input_file : str
        Path to the file that stores nucleotide details in (<nucleotide>, <matching symbol>) format at each position for each isomiR.
    output_file : str 
        Path to the file that stores the nucleotide or matching at each position for each isomiR.

    Example
    ----------- 
    ```
    The original file:                      (' ', ' '), (' ', ' '), ('g', '-'), ('a', '+'), ('u', '+'), ('c', '+'), ('c', '+'), ('u', '+'), ('g', '+')
    The output file (type = nt):            '','','g','a','u,'c','c','u','g' 
    The output file (type = templated):     '','','-','+','+,'+','+','+','+' 
    ```
    
    Returns 
    -------
    None. A new file that stores the nucleotide or matching at each position for each isomiR is generated. 
    """
    # Read input file
    templated_nt = pd.read_csv(input_file, low_memory=False)
    # Replace NA with ''
    templated_nt = templated_nt.fillna('')
    # Get position columns 
    cols = list(templated_nt.columns)
    del cols[0:4]
    # Loop over all rows 
    for index, r in templated_nt.iterrows():
        if r['is_pre'] == False:            
            for col in cols:
                vals = r[col]
                if vals != '':
                    first_val = ''
                    second_val = ''
                    if vals != "(' ', ' ')":
                        vals = vals.replace('(', '').replace(')', '').replace(' ', '').split(',')
                        first_val = vals[0]
                        second_val = vals[1]

                    if type == 'templated':
                        templated_nt.loc[index, col] = second_val
                    else: 
                        templated_nt.loc[index, col] = first_val  
    templated_nt.to_csv(output_file, index=False)

# Path to nt and templated alignment output files 
path_nt_templated_alignment_output_folder = sys.argv[1]
# Path to nt alignment output files 
path_nt_alignment_output_folder = sys.argv[2]
# Path to templated alignment output files 
path_templated_alignment_output_folder = sys.argv[3]
# List of group folders
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
        # Generate the templated alignment from nt templated file 
        split_nt_templated(f'{path_nt_templated_alignment_output_folder}/{group}/{rep_file}', f'{path_templated_alignment_output_folder}/{group}/{rep_file}', 'templated')