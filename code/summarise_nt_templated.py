import pandas as pd 
from collections import Counter
import os
import sys
from colorama import Fore, Style, init
init(autoreset=True)

def get_extension_values(r, max_nt_diff_5p, max_nt_diff_3p):
    """Get the matching symbol or nucleotide at extension positions (5'+1, 5'+2, 5'+3, ..., 3'+1, 3'+2, 3'+3,...) for each isomiR.

    Parameters
    ----------
    r : pandas.Series 
        A 1D array that stores details of an isomiRs such as sequence, matching symbol / nucleotides at each position. 
        Indexes correspond to column names, values correspond to cell values. 
    max_nt_diff_5p : int 
        The maximum number of nucleotide difference at 5' end across all isomiRs.
    max_nt_diff_3p : int 
        The maximum number of nucleotide difference at 3' end across all isomiRs.

    
    Example
    -----------
    ```
    r (nt alignment) : sja-bantam,   UGAGAUCGCGAUUAAAGCUGGU        ,False,,,,,u,g,a,g,a,u,c,g,c,g,a,u,u,a,a,a,g,c,u,g,g,u,,,,,,,,,,,
    max_nt_diff_5p : 3
    max_nt_diff_3p : 8
    Output (for nt alignment file) : ['a','g','u','a','a','g','c','u','g','g','u'] - the nucleotides at 5'+1,5'+2,5'+3,3'+1,3'+2,3'+3,3'+4,3'+5,3'+6,3'+7,3'+8.

    r (templated alignment) : sja-bantam,   UGGGAUCGCGAUUAAAGCUGGA        ,False,,,,,+,+,-,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,-,,,,,,,,,,,
    max_nt_diff_5p : 3
    max_nt_diff_3p : 8
    Output (for templated alignment file) : ['-','+','+','+','+','+','+','+','+','+','-'] - the matching symbols at 5'+1,5'+2,5'+3,3'+1,3'+2,3'+3,3'+4,3'+5,3'+6,3'+7,3'+8.
    ```
    
    Returns
    -------
    pandas.Series 
        A list of matching symbols or nucleotides at extension positions.
    """
    values = []
    pre_len = len(r['pre_seq'])
    for i in range(max_nt_diff_5p, 0, -1):
        values.append(r[str(i)])
    for i in range(pre_len - max_nt_diff_3p + 1, pre_len + 1):
        values.append(r[str(i)])
    return pd.Series(values)

def summarise_nt_alignment(path_nt_alignment_file, path_summarised_nt_alignment_file, max_nt_diff_5p, max_nt_diff_3p):
    """Calculate the nucleotide frequency at extension positions and save to the summarised nt alignment file.
    
    Parameters
    ----------
    path_nt_alignment_file : str 
        Path to nt alignment file. 
    path_summarised_nt_alignment_file : str
        Path to summarised nt alignment (for extension positions) file. 
    max_nt_diff_5p : int 
        The maximum number of nucleotide difference at 5' end across all isomiRs.
    max_nt_diff_3p : int 
        The maximum number of nucleotide difference at 3' end across all isomiRs.

    Example
    ----------- 
    ```
    Extension positions : 5'+1, 5'+2, 5'+3, 3'+1, 3'+2, 3'+3, 3'+4, 3'+5, 3'+6, 3'+7, 3'+8.
    sja-bantam :    u,g,g,u,g,u,a,c,g,c,a
    sja-bantam :    a,g,a,u,c,a,u,c,g,a,a
    sja-miR-10-3p : a,u,u,c,g,a,g,c,g,c,g
    sja-miR-10-3p : u,g,a,u,a,u,a,c,u,u,u


    Output :
    position | nucleotide | value
    5'+1     |          a |     2
    5'+1     |          u |     2        
    5'+1     |          c |     0
    5'+1     |          g |     0
    5'+2     |          a |     0
    5'+2     |          u |     1
    5'+2     |          c |     0
    5'+2     |          g |     3
    ....     |        ... |   ...
    ```

    Returns
    -------
    None. A summarised nt alignment file that stores the nucleotide frequency at extension positions is generated. 
    """
    # Read the nt alignment file 
    nt_alignment = pd.read_csv(path_nt_alignment_file)

    # Create a list of columns for extension positions at 5p
    extension_5p_cols = [ f"5'+{i + 1}" for i in range(max_nt_diff_5p)]
    # Create a list of columns for extension positions at 3p
    extension_3p_cols = [ f"3'+{i + 1}" for i in range(max_nt_diff_3p)]
    # Combine 5 extension cols and 3 extension cols
    extension_cols = extension_5p_cols + extension_3p_cols

    # Remove precursor rows 
    nt_alignment = nt_alignment[nt_alignment['is_pre'] == False]
    # For each isomiR/canonical, get nucleotide at 5e1, 5e2, 5e3,..., 3e1, 3e2,... 
    nt_alignment[extension_cols] = nt_alignment.apply(lambda r: get_extension_values(r, max_nt_diff_5p, max_nt_diff_3p), axis=1)
    # Create a dataframe that stores the nucleotide frequency at extension positions
    nt_summary = pd.DataFrame(columns=['position', 'nucleotide', 'value'])
    for col in extension_cols: 
        freq_counts = Counter(list(nt_alignment[col]))
        for nt in ['a', 'u', 'c', 'g']:
            nt_value = freq_counts[nt] if nt in freq_counts else 0
            nt_summary.loc[len(nt_summary.index)] = [col, nt, nt_value]
    nt_summary.to_csv(path_summarised_nt_alignment_file, index = False)

def summarise_templated_alignment(path_templated_alignment_file, path_summarised_templated_alignment_file, max_nt_diff_5p, max_nt_diff_3p):
    """Calculate the templated / nontemplated frequency at extension positions and save to the summarised templated alignment file.
    
    Parameters
    ----------
    path_templated_alignment_file : str 
        Path to templated alignment file. 
    path_summarised_templated_alignment_file : str
        Path to summarised templated alignment (for extension positions) file. 
    max_nt_diff_5p : int 
        The maximum number of nucleotide difference at 5' end across all isomiRs.
    max_nt_diff_3p : int 
        The maximum number of nucleotide difference at 3' end across all isomiRs.

    Example
    ------- 
    ```
    Extension positions : 5'+1, 5'+2, 5'+3, 3'+1, 3'+2, 3'+3, 3'+4, 3'+5, 3'+6, 3'+7, 3'+8.
    sja-bantam :     +,+,+,+,+,+,+,-,+,-,+
    sja-bantam :     -,-,+,+,+,+,+,+,+,-,-
    sja-miR-10-3p :  +,+,+,-,-,-,+,-,+,+,-
    sja-miR-10-3p :  -,+,+,+,+,+,+,+,-,-,-


    Output :
    position | templated   | value
    5'+1     |  Templated  |     2
    5'+1     |Nontemplated |     2        
    5'+2     |  Templated  |     3
    5'+2     |Nontemplated |     1
    5'+3     |  Templated  |     4
    5'+3     |Nontemplated |     0
    3'+1     |  Templated  |     3
    3'+2     |Nontemplated |     1
    ....     |        ...  |   ...
    ```

    Returns
    -------
    None. A summarised templated alignment file that stores the templated / nontemplated frequency at extension positions is generated. 
    """
    # Read the templated alignment file 
    templated_alignment = pd.read_csv(path_templated_alignment_file)
    # Create a list of columns for extension positions at 5p
    extension_5p_cols = [ f"5'+{i + 1}" for i in range(max_nt_diff_5p)]
    # Create a list of columns for extension positions at 3p
    extension_3p_cols = [ f"3'+{i + 1}" for i in range(max_nt_diff_3p)]
    # Combine 5 extension cols and 3 extension cols
    extension_cols = extension_5p_cols + extension_3p_cols

    # Remove precursor rows 
    templated_alignment = templated_alignment[templated_alignment['is_pre'] == False]
    # For each miRNA/isomiR, get values at 5e1, 5e2, 5e3,..., 3e1, 3e2,... 
    templated_alignment[extension_cols] = templated_alignment.apply(lambda r: get_extension_values(r, max_nt_diff_5p, max_nt_diff_3p), axis=1)
    # Create a dataframe that stores the templated / nontemplated frequency at extension positions
    templated_summary = pd.DataFrame(columns=['position', 'templated', 'value'])
    for col in extension_cols: 
        freq_counts = Counter(list(templated_alignment[col]))
        templated_value = freq_counts['+'] if '+' in freq_counts else 0
        untemplated_value = freq_counts['-'] if '-' in freq_counts else 0
        templated_summary.loc[len(templated_summary.index)] = [col, 'Templated', templated_value]
        templated_summary.loc[len(templated_summary.index)] = [col, 'Nontemplated', untemplated_value]
    templated_summary.to_csv(path_summarised_templated_alignment_file, index = False)

def summarise_templated_alignment_all(path_templated_alignment_file, path_summarised_templated_alignment_all_file, max_nt_diff_5p):
    """Calculate the templated / nontemplated frequency at all positions and save to the summarised templated alignment file. 

    Parameters
    ----------
    path_templated_alignment_file : str 
        Path to templated alignment file. 
    path_summarised_templated_alignment_all_file : str
        Path to summarised templated alignment for (all positions) file. 
    max_nt_diff_5p : str
        The maximum number of nucleotide difference at 5' end across all isomiRs.

    Example
    -------
    ```
    sja-bantam,   UGAGAUCGCGAUUAAAGCUGGUC       ,False,extended,,,,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,-,,,,,,,,,,
    sja-bantam,   UGAGAUCGCGAUUAAAGCUGGUAAU     ,False,extended,,,,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,-,-,+,,,,,,,,
    sja-miR-1,   UGGAAUGUGGCGAAGUAUGGUCA       ,False,extended,,,,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,-,,,,,,,,,,
    sja-miR-1,   UGAAAUGUGGCGAAGUAUGGUCU       ,False,extended,,,,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,+,,,,,,,,,,

    Output : 
    position | templated   | value
    5'+1     |  Templated  |     0
    5'+1     |Nontemplated |     0      
    5'+2     |  Templated  |     0
    5'+2     |Nontemplated |     0
    1        |  Templated  |     4
    1        |Nontemplated |     0
    ....     |     ...     |   ...
    23       |  Templated  |     1
    23       |Nontemplated |     3
    ....     |     ...     |   ...
    ```

    Returns 
    -------
    None. A summarised templated alignment file that stores the templated / nontemplated frequency at all positions is generated. 
    """
    # Read the templated alignment file 
    templated_alignment = pd.read_csv(path_templated_alignment_file)
    # Remove precursor rows 
    templated_alignment = templated_alignment[templated_alignment['is_pre'] == False]
    # Create a dataframe that stores the templated / nontemplated frequency at all positions
    templated_summary = pd.DataFrame(columns=['position', 'templated', 'value'])
    # Position columns 
    cols = list(templated_alignment.columns)
    del cols[0:4]

    # Loop through each position 
    for col in cols: 
        freq_counts = Counter(list(templated_alignment[col]))
        templated_value = freq_counts['+'] if '+' in freq_counts else 0
        untemplated_value = freq_counts['-'] if '-' in freq_counts else 0
        col = int(col)
        col = f"5'+{max_nt_diff_5p - col + 1}" if col <= max_nt_diff_5p else col - max_nt_diff_5p
        templated_summary.loc[len(templated_summary.index)] = [col, 'Templated', templated_value]
        templated_summary.loc[len(templated_summary.index)] = [col, 'Nontemplated', untemplated_value]
    templated_summary.to_csv(path_summarised_templated_alignment_all_file, index = False)
   
def run(
        path_nt_alignment_output_folder,
        path_templated_alignment_output_folder,
        path_summarised_nt_alignment_output_folder,
        path_summarised_templated_alignment_output_folder,
        path_summarised_templated_alignment_all_output_folder,
        path_precursors_output_folder):
    print(Fore.MAGENTA + "\nSummarising statistics for different types of variation ...")
    
    # List of group folders
    nt_group_folders = os.listdir(path_nt_alignment_output_folder)
    templated_group_folders = os.listdir(path_templated_alignment_output_folder)
    # Get precursor file 
    precursor_output_file = [file for file in os.listdir(path_precursors_output_folder) if '.csv' in file][0]
    # Get max nt difference at 5p 
    max_nt_diff_5p, max_nt_diff_3p = int(precursor_output_file.split('_')[0]), int(precursor_output_file.split('_')[1])

    for nt_group, templated_group in zip(nt_group_folders, templated_group_folders):
        # Get the list of nt alignment files of that group
        rep_nt_files = os.listdir(f'{path_nt_alignment_output_folder}/{nt_group}')
        # Get the list of templated alignment files of that group
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
                    
                                                                                                    