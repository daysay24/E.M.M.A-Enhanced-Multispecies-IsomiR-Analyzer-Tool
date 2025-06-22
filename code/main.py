import os 
import re
from colorama import Fore, Style, init
init(autoreset=True)

def print_menu():
    print(Fore.CYAN + "\nWelcome to the isomiR Analyzer Tool")
    print("=" * 40)
    print("Please select one of the following options:\n")
    print("  1. Analyse isomiRs of a species")
    print("  2. Visualise multi-species")
    print("  3. Stop")

def get_non_empty_str_value(msg): 
    while True:   
        print(msg, end='')  
        value = input()  
        if value: 
            return value 
        else: 
            print(Fore.RED + 'Must not be empty! Please try again.')
        
def get_yes_no_value(msg):
    while True: 
        print(msg, end='')
        value = input()
        if value == 'Y' or value == 'y':
            return True
        elif value == 'N' or value == 'n':
            return False
        else: 
            print(Fore.RED + 'Must be Y/y or N/n! Please try again.')

def get_numeric_value(msg):
    while True: 
        print(msg, end='')
        value = input() 

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
        
def analyse_isomirs():
    print(Fore.CYAN + "\nStep 1: Input Species Information")
    print("-" * 35)

    root_folder = os.getcwd().replace('code', '')

    input_root_folder = root_folder + 'input'
    print("Input root folder is:", Fore.GREEN + input_root_folder)

    species_code = get_non_empty_str_value('Species code: ')

    species_name = get_non_empty_str_value('Species name (for labeling in visualisation): ')

    print(Fore.CYAN + "\nStep 2: Checking Required Files")
    print("-" * 35)
    print(f"Input folder for species '{species_name} ({species_code})' is located at:")
    print(Fore.GREEN + f"  {input_root_folder}/{species_code}")
    print("\nPlease ensure this folder contains the following required files:")
    print("  - genome.fa")
    print("  - UTR.fa")
    print("  - miRNA_annotation.gff3")
    print("  - isomiR_SEA_outputs/")

    print(Fore.CYAN + "\nStep 3: Parameters and Output")
    print("-" * 35)

    read_count_thres = get_numeric_value('Read count threshold (Press Enter to use default 10): ')
    if not read_count_thres:
        read_count_thres = 10

    is_mirbase_gff = get_yes_no_value('Is miRNA_annotation.gff3 from mirbase (Y/N) ?:')
    
    is_match_chr_names = get_yes_no_value('Is match gff to genome required (Y/N) ?:')

    output_root_folder = root_folder + 'output'
    print("Output root folder is:", Fore.GREEN + f"{output_root_folder}/{species_code}")

def visualise_isomirs(): 
    print(Fore.MAGENTA + "\n[Visualising multi-species ... work in progress]\n")

def stop(): 
    print(Fore.GREEN + "\nGoodbye!\n")

while True: 
    print_menu()
    choice = input(Fore.YELLOW + "\nEnter your choice: ").strip()

    if choice == '1':
        analyse_isomirs()
    elif choice == '2':
        visualise_isomirs()
    elif choice == '3': 
        stop()
        break
    else: 
        print(Fore.RED + "\nInvalid choice. Please enter 1, 2, or 3.")