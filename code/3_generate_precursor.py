import pandas as pd 
import subprocess
import re
import os
import re
import requests 
import sys

def extract_miR_info(row):
    details = row['details'].split(';')
    if (row['type'] == "miRNA_primary_transcript"):
        return pd.Series([details[1].split('=')[1], details[2].split('=')[1], ''])
    else:
        return pd.Series([details[1].split('=')[1], details[2].split('=')[1], details[3].split('=')[1]])
    
def create_bed_coordinates(is_built_in_genome, genomic_coordinates, max_nt_diff_5p, max_nt_diff_3p, path_precursor_output_folder):
    shift = 0 if is_built_in_genome == True else 1
    # change the canoical to lower (from the position 0-max5p - 1 and end-max3p-1)
    with open(f'{path_precursor_output_folder}/extended_precursor_coords.bed', "w+") as bed_file: 
        for _, r in genomic_coordinates.iterrows(): 
            if r['strand'] == '+':
                bed_file.write(f"{r['chr']}\t{int(r['start']) - max_nt_diff_5p - shift}\t{int(r['end']) + max_nt_diff_3p}\t{r['name']}\t1\t{r['strand']}\n")
            else: 
                bed_file.write(f"{r['chr']}\t{int(r['start']) - max_nt_diff_3p - shift}\t{int(r['end']) + max_nt_diff_5p}\t{r['name']}\t1\t{r['strand']}\n")

def get_extended_miRNA_coordinates(is_mirbase_gff, is_built_in_genome, max_nt_diff_5p, max_nt_diff_3p, path_precursors_output_folder, path_genomic_file, path_coords_file, species=None):
    # read the annotation file 
    if is_mirbase_gff == True: 
        genomic_coordinates = pd.read_csv(path_coords_file, sep='\t', header=None, skiprows=[0,1,2,3,4,5,6,7,8,9,10,11,12], names=["chr", "unknown1", "type", "start", "end", "unknown2", "strand", "unknown3", "details"])
        genomic_coordinates[['id', 'name', 'derives_from']] = genomic_coordinates.apply(lambda r: extract_miR_info(r), axis = 1)
        genomic_coordinates = genomic_coordinates[(genomic_coordinates['type'] == 'miRNA')]
        genomic_coordinates = genomic_coordinates[['chr', 'name', 'start', 'end', 'strand']]
    else:
        genomic_coordinates = pd.read_excel(path_coords_file)

    # create bed coords 
    create_bed_coordinates(is_built_in_genome, genomic_coordinates, max_nt_diff_5p, max_nt_diff_3p, path_precursors_output_folder)   

    if is_built_in_genome:
        # miRNAs 
        miRNAs = []
        # extended precursor sequences 
        extended_precursor_seqs = []
        # Read the bed file 
        extended_precursor_ranges = pd.read_csv(f'{path_precursors_output_folder}/extended_precursor_coords.bed', sep='\t', names=['chr', 'start', 'end', 'mir_name', 'no', 'strand'])

        # Loop through each miRNA and call the REST API to get the extended sequence 
        for i, row in extended_precursor_ranges.iterrows():
            base_query = f'http://togows.org/api/ucsc/{species}/{row["chr"]}:'           
            query = f'{base_query}{row["start"]}-{row["end"]}' if row['strand'] == '+' else f'{base_query}{row["end"]}-{row["start"]}'
            r = requests.get(query)
            extended_precursor_seq = r.text.upper().replace("T", "U")
            miRNAs.append(row['mir_name'])
            extended_precursor_seqs.append(extended_precursor_seq)
        pd.DataFrame({'mir_name': miRNAs, 'extended_precursor_seq': extended_precursor_seqs}).to_csv(f'{path_precursors_output_folder}/{max_nt_diff_5p}_{max_nt_diff_3p}_extended_precursor_seqs.csv', index=False)
    else:
        # run BEDTools to get extended precursor from genomic coordinates
        abs_path_to_genome_file = os.path.abspath(path_genomic_file)
        abs_path_to_extended_coords_bed_file = os.path.abspath(f'{path_precursors_output_folder}/extended_precursor_coords.bed') 
        abs_path_to_extracted_precursor_seqs_file = os.path.abspath(f'{path_precursors_output_folder}/extended_precursor_seqs.fa')
        cmd = f'bedtools getfasta -fi {abs_path_to_genome_file} -bed {abs_path_to_extended_coords_bed_file} -fo {abs_path_to_extracted_precursor_seqs_file} -s -nameOnly'
        subprocess.call(cmd, shell=True)
        # From the output file of BEDTools, save all miRNA names and their extended precursor sequences in a csv file 
        mir_names = []
        extended_precursor_seqs = []

        for line in open(abs_path_to_extracted_precursor_seqs_file).readlines():
            if '>' in line:
                mir_names.append(re.search(r'>(.*?)\(', line).group(1))
            else:
                extended_precursor_seq = line.replace('T', 'U').strip()
                miRNA_seq = extended_precursor_seq[max_nt_diff_5p:len(line) - 1 - max_nt_diff_3p]
                extended_precursor_seq = extended_precursor_seq.replace(miRNA_seq, miRNA_seq.lower())
                extended_precursor_seqs.append(extended_precursor_seq)
        pd.DataFrame({'mir_name': mir_names, 'extended_precursor_seq': extended_precursor_seqs}).to_csv(f'{path_precursors_output_folder}/{max_nt_diff_5p}_{max_nt_diff_3p}_extended_precursor_seqs.csv', index=False)

# Path to the summarised outputs folder 
path_summarised_output_folder = sys.argv[1]
# Path to the precursor folder 
path_precursors_output_folder = sys.argv[2]
# Path to the genomic reference file 
path_genomic_file = sys.argv[3]
# Path to the miRNA coordinates file 
path_coords_file = sys.argv[4]
# Species code 
species = sys.argv[5]
# Is user using mirbase gff file 
is_mirbase_gff = True if sys.argv[6] == 'True' else False

print(path_summarised_output_folder, path_precursors_output_folder, path_genomic_file, path_coords_file, species, is_mirbase_gff)

# Is user using built_in genome 
is_built_in_genome = True if sys.argv[7] == 'True' else False


# Max nt difference at 5p and 3p 
max_nt_diff_5p, max_nt_diff_3p = 0, 0
# create folder if not exists
if not os.path.exists(path_precursors_output_folder):
    os.makedirs(path_precursors_output_folder)

# Loop through each group
group_folders = os.listdir(path_summarised_output_folder)
for group in group_folders:
    # Get the list of replicate files
    rep_files = os.listdir(f'{path_summarised_output_folder}/{group}')
    # Loop through each replicate of that group 
    for rep_file in rep_files:
        # Get the replicate name 
        rep_name = rep_file.split('.')[0]
        # Read summarised isomiRs file of that replicate 
        summarised_isomiRs = pd.read_csv(f'{path_summarised_output_folder}/{group}/{rep_file}', encoding='latin-1')
        # Update max nt difference at 5p and 3p if necessary
        if max(summarised_isomiRs['5p_nt_diff']) > max_nt_diff_5p: 
            max_nt_diff_5p = max(summarised_isomiRs['5p_nt_diff'])
        if max(summarised_isomiRs['3p_nt_diff']) > max_nt_diff_3p:
            max_nt_diff_3p = max(summarised_isomiRs['3p_nt_diff'])
get_extended_miRNA_coordinates(is_mirbase_gff, is_built_in_genome, max_nt_diff_5p, max_nt_diff_3p, path_precursors_output_folder, path_genomic_file, path_coords_file, species)
