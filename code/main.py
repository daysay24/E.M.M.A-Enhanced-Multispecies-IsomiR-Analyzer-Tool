import os 
import time
import subprocess
import re
import signal
import pandas as pd
import summarise_isomir_sea 
import avg_summarised_isomirs 
import generate_precursor
import nt_templated
import split_nt_templated
import summarise_nt_templated
import avg_summarised_nt_templated
import process_graph_data
from colorama import Fore, Style, init
init(autoreset=True)

def print_menu():
    print(Fore.CYAN + "\nWelcome to the isomiR Analyzer Tool")
    print("=" * 40)
    print("Please select one of the following options:\n")
    print("  1. Analyse isomiRs of a species")
    print("  2. Visualise multi-species")
    print("  3. Stop")
    print("Or press Ctrl+C to force quit.")

def get_non_empty_str_value(msg): 
    while True:   
        value = input(Fore.YELLOW + msg).strip()  
        if value: 
            return value 
        else: 
            print(Fore.RED + 'Must not be empty! Please try again.')
        
def get_yes_no_value(msg):
    while True: 
        value = input(Fore.YELLOW + msg).strip()
        if value == 'Y' or value == 'y':
            return True
        elif value == 'N' or value == 'n':
            return False
        else: 
            print(Fore.RED + 'Must be Y/y or N/n! Please try again.')

def get_numeric_value(msg):
    while True: 
        value = input(Fore.YELLOW + msg).strip()

        if not value: 
            return None

        if re.match(r'^[+-]?[0-9]+$', value): 
            value = int(value)
            if value < 0: 
                print(Fore.RED + 'Must be a positive value! Please try again.')
            else: 
                return value
        else: 
            print(Fore.RED + 'Must be an integer! Please try again.')

def check_input_files_exist(input_folder):
    if not os.path.exists(f"{input_folder}"):
        raise FileNotFoundError(f"{input_folder}/ not found! Please add the required folder and try again.")
    if not os.path.exists(f"{input_folder}/genomic.fa"):
        raise FileNotFoundError(f"{input_folder}/genomic.fa not found! Please add the required file and try again.")
    if not os.path.exists(f"{input_folder}/UTR.fa"):
        raise FileNotFoundError(f"{input_folder}/UTR.fa not found! Please add the required file and try again.")
    if not os.path.exists(f"{input_folder}/miRNA_annotation.gff3") and not os.path.exists(f"{input_folder}/miRNA_annotation.xlsx"):
        raise FileNotFoundError(f"Both {input_folder}/miRNA_annotation.gff3 and {input_folder}/miRNA_annotation.xlsx not found! Please add either required file and try again.")
    if not os.path.exists(f"{input_folder}/isomiR-SEA_outputs/"):
        raise FileNotFoundError(f"{input_folder}/isomiR-SEA_outputs/ not found! Please add the required folder and try again.")

def update_metadata_file(species_code, species_name, input_folder, root_folder):
    group_folders = os.listdir(input_folder + '/isomiR-SEA_outputs/')
    n_groups = len(group_folders)
    metadata_path = root_folder + '/dashboard/metadata.csv'
    species_medadata_df = pd.DataFrame({
        'species': [species_code] * n_groups, 
        'alias': [species_name] * n_groups,
        'group': group_folders
    })
    if not os.path.exists(metadata_path):
        species_medadata_df.to_csv(metadata_path, index = False)
    else: 
        current_medadata_df = pd.read_csv(metadata_path)
        current_medadata_df = current_medadata_df[(current_medadata_df['species'] != species_code) & (current_medadata_df['alias'] != species_name)]
        all_medadata_df = pd.concat([current_medadata_df, species_medadata_df], ignore_index=True)
        all_medadata_df.to_csv(metadata_path, index = False)
    
def get_analyse_isomirs_info(): 
    print(Fore.CYAN + "\nStep 1: Input Species Information")
    print("-" * 35)

    root_folder = os.getcwd().replace('code', '')

    species_code = get_non_empty_str_value('Species code: ')

    species_name = get_non_empty_str_value('Species name (for labeling in visualisation): ')

    input_root_folder = root_folder + 'input'
    input_folder = f"{input_root_folder}/{species_code}"
    print("Input folder is:", Fore.GREEN + input_folder)

    print(Fore.CYAN + "\nStep 2: Checking Required Files")
    print("-" * 35)
    print(f"Input folder for species '{species_name} ({species_code})' is located at:")
    print(Fore.GREEN + f"  {input_root_folder}/{species_code}")
    print("\nPlease ensure this folder contains the following required files:")
    print("  - genomic.fa")
    print("  - UTR.fa")
    print("  - miRNA_annotation.gff3 or miRNA_annotation.xlsx")
    print("  - isomiR-SEA_outputs/")

    print(Fore.CYAN + "\nStep 3: Parameters and Output")
    print("-" * 35)

    read_count_thres = get_numeric_value('Read count threshold (Press Enter to use default 10): ')
    if not read_count_thres:
        read_count_thres = 10

    is_mirbase_gff = get_yes_no_value('Is miRNA_annotation file from mirbase (Y/N) ?:')
    
    is_match_chr_names = get_yes_no_value('Is match gff to genome required (Y/N) ?:')

    output_root_folder = root_folder + 'output'
    output_folder = f"{output_root_folder}/{species_code}"
    print("Output folder is:", Fore.GREEN + output_folder)

    return root_folder, input_folder, species_code, species_name, read_count_thres, is_mirbase_gff, is_match_chr_names, output_folder

def analyse_isomirs():
    root_folder, input_folder, species_code, species_name, read_count_thres, is_mirbase_gff, is_match_chr_names, output_folder = get_analyse_isomirs_info()

    path_genomic_file = input_folder + '/genomic.fa'
    path_coords_file = input_folder + '/miRNA_annotation.gff3' if is_mirbase_gff else input_folder + '/miRNA_annotation.xlsx'
    path_raw_output_folder = input_folder + '/isomiR-SEA_outputs'
    path_summarised_output_folder = output_folder + '/1_summarised_isomiRs'
    path_avg_replicate_output_folder = output_folder + '/2_avg_replicate_isomiRs'
    path_precursors_output_folder = output_folder + '/3_precursors'
    path_nt_templated_alignment_output_folder = output_folder + '/4_nt_templated_alignment'
    path_nt_alignment_output_folder = output_folder + '/5_nt_alignment'
    path_templated_alignment_output_folder = output_folder + '/5_templated_alignment'
    path_summarised_nt_alignment_output_folder = output_folder + '/6_summarised_nt_alignment'
    path_summarised_templated_alignment_output_folder = output_folder + '/6_summarised_templated_alignment'
    path_summarised_templated_alignment_all_output_folder = output_folder + '/6_summarised_templated_alignment_all'
    path_avg_summarised_nt_alignment_output_folder = output_folder + '/7_avg_summarised_nt_alignment'
    path_avg_summarised_templated_alignment_output_folder = output_folder + '/7_avg_summarised_templated_alignment'
    path_avg_summarised_templated_alignment_all_output_folder = output_folder + '/7_avg_summarised_templated_alignment_all'
    path_graph_processed_data_folder = output_folder + '/8_graph_processed_data/'

    try: 
        check_input_files_exist(input_folder)
        summarise_isomir_sea.run(
            path_raw_output_folder, 
            path_summarised_output_folder, 
            read_count_thres)
        avg_summarised_isomirs.run(
            path_summarised_output_folder, 
            path_avg_replicate_output_folder)
        generate_precursor.run(
            path_summarised_output_folder, 
            path_precursors_output_folder, 
            path_genomic_file, 
            path_coords_file, 
            is_mirbase_gff, 
            is_match_chr_names)
        nt_templated.run(
            path_summarised_output_folder, 
            path_precursors_output_folder, 
            path_nt_templated_alignment_output_folder)
        split_nt_templated.run(
            path_nt_templated_alignment_output_folder, 
            path_nt_alignment_output_folder, 
            path_templated_alignment_output_folder)
        summarise_nt_templated.run(
            path_nt_alignment_output_folder,
            path_templated_alignment_output_folder,
            path_summarised_nt_alignment_output_folder,
            path_summarised_templated_alignment_output_folder,
            path_summarised_templated_alignment_all_output_folder,
            path_precursors_output_folder)
        avg_summarised_nt_templated.run(
            path_summarised_nt_alignment_output_folder,
            path_summarised_templated_alignment_output_folder,
            path_summarised_templated_alignment_all_output_folder,
            path_avg_summarised_nt_alignment_output_folder,
            path_avg_summarised_templated_alignment_output_folder,
            path_avg_summarised_templated_alignment_all_output_folder)
        process_graph_data.run(
            path_avg_replicate_output_folder,
            path_avg_summarised_templated_alignment_output_folder,
            path_avg_summarised_nt_alignment_output_folder,
            path_avg_summarised_templated_alignment_all_output_folder,
            path_graph_processed_data_folder)
        update_metadata_file(species_code, species_name, input_folder, root_folder)
    except Exception as e: 
        print(f'Analyse isomiRs of {species_name} ({species_code}) failed due to: {e}')

def visualise_isomirs(): 
    print(Fore.MAGENTA + "\nVisualising multi-species work in progress ...\n")

    with open("../dashboard/error.log", "w") as log_file:
        process = subprocess.Popen(
            ["python", "../dashboard/app.py"],
            stdout=log_file,
            stderr=log_file,
            # preexec_fn=os.setpgrp
        )
    time.sleep(5)
    return_code = process.poll()
    if return_code == None: 
        print(Fore.MAGENTA + "\nOpen http://127.0.0.1:8050/ in your browser\n")
        return process    
    else:
        print(Fore.RED + "\nSome thing went wrong. Please check the /dashboard/error.log file !")

def kill_process(process):
    os.killpg(os.getpgid(process.pid), signal.SIGKILL)

def is_visualisation_process_stopped(process): 
    is_stopping_visualisation = get_yes_no_value(Fore.RED + "\nYou can't add more species while the visualisation is running. Do you want to stop it first ?: ")
    if is_stopping_visualisation: 
        kill_process(process)
        return True
    else:
        return False

if __name__ == "__main__":
    process = None
    while True: 
        print_menu()
        choice = input(Fore.YELLOW + "\nEnter your choice: ").strip()

        if choice == '1':
            if not process: 
                analyse_isomirs()
            else: 
                if is_visualisation_process_stopped(process):
                    process = None
                    analyse_isomirs()
        elif choice == '2':
            if process: 
                print(Fore.MAGENTA + "\nVisualising multi-species work in progress\n")
            else:
                process = visualise_isomirs()
        elif choice == '3':
            if process: 
                kill_process(process)
            print(Fore.GREEN + "\nExiting the program...\n")
            break
        else: 
            print(Fore.RED + "\nInvalid choice. Please enter 1, 2, or 3.")