#!/bin/bash

# Delete all .DS_store files 
find . -name ".DS_Store" -type f -delete

# Create an output folder if not exist. If already exist, delete all files 
path_output_folder='./output'
# Check if the output directory exists
if [ -d "$path_output_folder" ]; then
  # If it exists, delete all files inside it
  rm -rf "${path_output_folder:?}"/*
else
  # If it does not exist, create the directory
  mkdir "$path_output_folder"
fi

# Paths 
path_genomic_file=./test/s.jap/s.jap.fa # change this if necessary
path_coords_file=./test/s.jap/s.jap.xlsx # change this if necessary
path_raw_output_folder=./test/s.jap/isomiR-SEA_outputs # change this if necessary
path_summarised_output_folder=./output/1_summarised_isomiRs
path_avg_replicate_output_folder=./output/2_avg_replicate_isomiRs
path_precursors_output_folder=./output/3_precursors
path_nt_templated_alignment_output_folder=./output/4_nt_templated_alignment
path_nt_alignment_output_folder=./output/5_nt_alignment
path_templated_alignment_output_folder=./output/5_templated_alignment
path_summarised_nt_alignment_output_folder=./output/6_summarised_nt_alignment
path_summarised_templated_alignment_output_folder=./output/6_summarised_templated_alignment
path_summarised_templated_alignment_all_output_folder=./output//6_summarised_templated_alignment_all
path_avg_summarised_nt_alignment_output_folder=./output//7_avg_summarised_nt_alignment
path_avg_summarised_templated_alignment_output_folder=./output/7_avg_summarised_templated_alignment
path_avg_summarised_templated_alignment_all_output_folder=./output/7_avg_summarised_templated_alignment_all
path_graph_processed_data_folder=./output/8_graph_processed_data/
path_graphs_folder=./output/graphs/

# Summarise isomiRs 
read_count_threshold=10
python ./code/1_summarise_isomiR_SEA.py $path_raw_output_folder $path_summarised_output_folder $read_count_threshold

# Average rep isomiRs 
python ./code/2_avg_summarised_isomiRs.py $path_summarised_output_folder $path_avg_replicate_output_folder

# Get precursors 
is_mirbase_gff=False # change this if necessary
match_chr_names=True # change this if necessary
python ./code/3_generate_precursor.py $path_summarised_output_folder $path_precursors_output_folder $path_genomic_file $path_coords_file $is_mirbase_gff $match_chr_names

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
Rscript ./code/9_visualise.R $path_graphs_folder $path_graph_processed_data_folder
