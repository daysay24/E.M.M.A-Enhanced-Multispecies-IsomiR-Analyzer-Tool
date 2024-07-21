#!/bin/bash

# Paths 
path_genomic_data_folder=./data/0_genomic_data
path_raw_output_folder=./data/0_isomiR-SEA_isomiRs
path_summarised_output_folder=./data/1_summarised_isomiRs
path_avg_replicate_output_folder=./data/2_avg_replicate_isomiRs
path_precursors_output_folder=./data/3_precursors
path_nt_templated_alignment_output_folder=./data/4_nt_templated_alignment
path_nt_alignment_output_folder=./data/5_nt_alignment
path_templated_alignment_output_folder=./data/5_templated_alignment
path_summarised_nt_alignment_output_folder=./data/6_summarised_nt_alignment
path_summarised_templated_alignment_output_folder=./data/6_summarised_templated_alignment
path_summarised_templated_alignment_all_output_folder=./data/6_summarised_templated_alignment_all
path_avg_summarised_nt_alignment_output_folder=./data/7_avg_summarised_nt_alignment
path_avg_summarised_templated_alignment_output_folder=./data/7_avg_summarised_templated_alignment
path_avg_summarised_templated_alignment_all_output_folder=./data/7_avg_summarised_templated_alignment_all
path_graph_processed_data_folder=./data/8_graph_processed_data/
path_graphs=./data/graphs/

# Summarise isomiRs 
python ./code/1_summarise_isomiR_SEA.py $path_raw_output_folder $path_summarised_output_folder

# Average rep isomiRs 
python ./code/2_avg_summarised_isomiRs.py $path_summarised_output_folder $path_avg_replicate_output_folder

# Get precursors 
# species=mm10
# is_mirbase_gff=True
# is_built_in_genome=True
# python ./code/3_generate_precursor $path_summarised_output_folder $path_precursors_output_folder $path_genomic_data_folder $species $is_mirbase_gff $is_built_in_genome

# Align nt and templated 
python ./code/4_nt_templated.py $path_summarised_output_folder $path_precursors_output_folder $path_nt_templated_alignment_output_folder

# Split nt and templated alignment
python ./code/5_split_nt_templated.py $path_nt_templated_alignment_output_folder $path_nt_alignment_output_folder $path_templated_alignment_output_folder

# Summarise nt and templated 
python ./code/6_summarise_nt_templated.py $path_nt_alignment_output_folder $path_templated_alignment_output_folder $path_summarised_nt_alignment_output_folder $path_summarised_templated_alignment_output_folder $path_summarised_templated_alignment_all_output_folder $path_precursors_output_folder

# Average summarised nt and templated
python ./code/7_avg_nt_templated.py $path_summarised_nt_alignment_output_folder $path_summarised_templated_alignment_output_folder $path_summarised_templated_alignment_all_output_folder $path_avg_summarised_nt_alignment_output_folder $path_avg_summarised_templated_alignment_output_folder $path_avg_summarised_templated_alignment_all_output_folder

# Process graph data 
python ./code/8_process_graph_data.py $path_avg_replicate_output_folder $path_avg_summarised_templated_alignment_output_folder $path_avg_summarised_nt_alignment_output_folder $path_avg_summarised_templated_alignment_all_output_folder $path_graph_processed_data_folder

# Generate graphs 
Rscript ./code/9_visualise.R $path_graphs $path_graph_processed_data_folder
