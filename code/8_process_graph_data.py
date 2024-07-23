import pandas as pd 
import os
import sys
# Graph 1: summarise data for showing miRNAs and isomiRs total reads (rpm) and relative abundance as percentage of total reads, across stages 
def process_graph1_data():
    # All groups df
    all_group_df = pd.DataFrame()
    # Path to the average isomiR files 
    path_avg_replicate_output_folder = sys.argv[1]
    # Path to graph processed data files 
    path_graph_processed_data_folder = sys.argv[5]
    # Loop through averaged summarised isomiRs of each group
    avg_files = os.listdir(path_avg_replicate_output_folder)
    for avg_file in avg_files:
        # Read averaged summarised isomiRs of each group
        avg_summarised_isomiRs_df = pd.read_csv(f'{path_avg_replicate_output_folder}/{avg_file}')
        # Select a subset of important columms 
        type_df = avg_summarised_isomiRs_df[['grouped_type', 'rpm']]
        # Add a new column type: if grouped_type is canonical, type is miR. otherwise, type is miR
        type_df['type'] = type_df['grouped_type'].apply(lambda gt: 'Canonical' if gt == 'Canonical' else 'IsomiR')
        # Drop grouped_type column 
        type_df = type_df.drop(columns=['grouped_type'])
        # Group by type column and sum the rpm values
        type_df = type_df.groupby(by='type').sum().reset_index()
        # Add relative abundance column 
        sum_rpm = type_df['rpm'].sum()
        type_df['relative_abundance'] = type_df['rpm'] / sum_rpm * 100
        # Add group column 
        type_df['group'] = avg_file.split('.')[0]
        
        if all_group_df.empty: 
            all_group_df = type_df
        else:
            all_group_df = pd.concat([all_group_df, type_df], ignore_index=True)
    all_group_df.to_csv(f'{path_graph_processed_data_folder}/graph_1_data.csv', index=False)

# Graph 2: summarise data for showing the relative abundance as percentage of total reads of types (5p, 3p, both, canonical, others) across stages
def process_graph2_data():
    # All groups df
    all_group_df = pd.DataFrame()
    # Path to the average isomiR files 
    path_avg_replicate_output_folder = sys.argv[1]
    # Path to graph processed data files 
    path_graph_processed_data_folder = sys.argv[5]
    # Loop through averaged summarised isomiRs of each group
    avg_files = os.listdir(path_avg_replicate_output_folder)
    for avg_file in avg_files:
        # Read averaged summarised isomiRs of each group
        avg_summarised_isomiRs_df = pd.read_csv(f'{path_avg_replicate_output_folder}/{avg_file}')
        # Select a subset of important columms 
        grouped_type_df = avg_summarised_isomiRs_df[['grouped_type', 'rpm', 'unique_tag']]
        # Group by grouped_type column and sum the rpm values and count unique tags 
        grouped_type_df = grouped_type_df.groupby('grouped_type').agg(rpm=('rpm', 'sum'), unique_tag=('unique_tag', 'sum')).reset_index()
        # Add group column 
        grouped_type_df['group'] = avg_file.split('.')[0]
        
        if all_group_df.empty: 
            all_group_df = grouped_type_df
        else:
            all_group_df = pd.concat([all_group_df, grouped_type_df], ignore_index=True)
    all_group_df.to_csv(f'{path_graph_processed_data_folder}/graph_2_data.csv', index=False)

# summarise data for showing proportions of 5p/3p addition/truncation at different positions (3e1, 3e2, 3e3,..., 3t1, 3t2, 5e1, 5e2,..., 5t1, 5t2, 5t3, ...) across stages 
def process_graph3_data():
    # All groups df
    all_group_df = pd.DataFrame()
    # Path to the average isomiR files 
    path_avg_replicate_output_folder = sys.argv[1]
    # Path to graph processed data files 
    path_graph_processed_data_folder = sys.argv[5]
    # Loop through averaged summarised isomiRs of each group
    avg_files = os.listdir(path_avg_replicate_output_folder)
    for avg_file in avg_files:
        # Read averaged summarised isomiRs of each group
        avg_summarised_isomiRs_df = pd.read_csv(f'{path_avg_replicate_output_folder}/{avg_file}')
        # Select a subset of important columms 
        type_nt_df = avg_summarised_isomiRs_df[['type_nt', 'rpm', 'unique_tag', 'grouped_type']]
        # Select isomiR 3p or 5p 
        type_nt_df = type_nt_df[type_nt_df['grouped_type'].isin(["5'isomiR", "3'isomiR"])]
        # Group by type_nt column and sum the rpm values and count unique tags
        type_nt_df = type_nt_df.groupby(['type_nt', 'grouped_type']).agg(rpm=('rpm', 'sum'), unique_tag=('unique_tag', 'sum')).reset_index()
        # Add group column 
        type_nt_df['group'] = avg_file.split('.')[0]
        
        if all_group_df.empty: 
            all_group_df = type_nt_df
        else:
            all_group_df = pd.concat([all_group_df, type_nt_df], ignore_index=True)
    all_group_df.to_csv(f'{path_graph_processed_data_folder}/graph_3_data.csv', index=False)

# summarise data for showing proportion of templated vs untemplated at addition positions in different groups. 
def process_graph4_data():
    # All groups df
    all_group_df = pd.DataFrame()
    # Path to the average templated summarised alignment files 
    path_avg_summarised_templated_alignment_output_folder = sys.argv[2]
    # Path to graph processed data files 
    path_graph_processed_data_folder = sys.argv[5]
    # Loop through averaged templated summarised alignment file of each group
    avg_files = os.listdir(path_avg_summarised_templated_alignment_output_folder)
    for avg_file in avg_files:
        # Read averaged templated summarised alignment file of each group
        avg_templated_summarised_alignment_df = pd.read_csv(f'{path_avg_summarised_templated_alignment_output_folder}/{avg_file}', dtype={'position': 'str'})
        # Add group column 
        avg_templated_summarised_alignment_df['group'] = avg_file.split('.')[0]
        
        if all_group_df.empty: 
            all_group_df = avg_templated_summarised_alignment_df
        else:
            all_group_df = pd.concat([all_group_df, avg_templated_summarised_alignment_df], ignore_index=True)
    all_group_df.to_csv(f'{path_graph_processed_data_folder}/graph_4_data.csv', index=False)

# summarise data for showing proportion of nucleotides (A, U, C, G) at addition positions in different groups.
def process_graph5_data():
    # All groups df
    all_group_df = pd.DataFrame()
    # Path to the average nt summarised alignment files 
    path_avg_summarised_nt_alignment_output_folder = sys.argv[3]
    # Path to graph processed data files 
    path_graph_processed_data_folder = sys.argv[5]
    # Loop through averaged nt summarised alignment file of each group
    avg_files = os.listdir(path_avg_summarised_nt_alignment_output_folder)
    for avg_file in avg_files:
        # Read averaged nt summarised alignment file of each group
        avg_nt_summarised_alignment_df = pd.read_csv(f'{path_avg_summarised_nt_alignment_output_folder}/{avg_file}',  dtype={'position': 'str'})
        # Add group column 
        avg_nt_summarised_alignment_df['group'] = avg_file.split('.')[0]
        
        if all_group_df.empty: 
            all_group_df = avg_nt_summarised_alignment_df
        else:
            all_group_df = pd.concat([all_group_df, avg_nt_summarised_alignment_df], ignore_index=True)
    all_group_df.to_csv(f'{path_graph_processed_data_folder}/graph_5_data.csv', index=False)

def process_graph6_data():
    # All groups df
    all_group_df = pd.DataFrame()
    # Path to the average templated summarised alignment files 
    path_avg_summarised_templated_alignment_all_output_folder = sys.argv[4]
    # Path to graph processed data files 
    path_graph_processed_data_folder = sys.argv[5]
    # Loop through averaged templated summarised alignment all file of each group
    avg_files = os.listdir(path_avg_summarised_templated_alignment_all_output_folder)
    for avg_file in avg_files:
        # Read averaged templated summarised alignment file of each group
        avg_templated_summarised_alignment_all_df = pd.read_csv(f'{path_avg_summarised_templated_alignment_all_output_folder}/{avg_file}',  dtype={'position': 'str'})
        # Add group column 
        avg_templated_summarised_alignment_all_df['group'] = avg_file.split('.')[0]
        
        if all_group_df.empty: 
            all_group_df = avg_templated_summarised_alignment_all_df
        else:
            all_group_df = pd.concat([all_group_df, avg_templated_summarised_alignment_all_df], ignore_index=True)
    all_group_df.to_csv(f'{path_graph_processed_data_folder}/graph_6_data.csv', index=False)


if not os.path.exists(sys.argv[5]):
    os.makedirs(sys.argv[5])

process_graph1_data()
process_graph2_data()
process_graph3_data()
process_graph4_data()
process_graph5_data()
process_graph6_data()