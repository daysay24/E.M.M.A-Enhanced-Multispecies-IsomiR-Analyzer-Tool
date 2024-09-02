import pandas as pd 
import subprocess
import re
import os
import re
import sys
import fileinput


single_chr_species = ['hhv507799', 'hehcmvcg', 'ksu75698']

def extract_miR_info(row):
    details = row['details'].split(';')
    if (row['type'] == "miRNA_primary_transcript"):
        return pd.Series([details[1].split('=')[1], details[2].split('=')[1], ''])
    else:
        return pd.Series([details[1].split('=')[1], details[2].split('=')[1], details[3].split('=')[1]])
    
def create_bed_coordinates(genomic_coordinates, max_nt_diff_5p, max_nt_diff_3p, path_precursor_output_folder):
    shift = 1
    # change the canoical to lower (from the position 0-max5p - 1 and end-max3p-1)
    with open(f'{path_precursor_output_folder}/extended_precursor_coords.bed', "w+") as bed_file: 
        for _, r in genomic_coordinates.iterrows(): 
            if r['strand'] == '+':
                bed_file.write(f"{r['chr']}\t{int(r['start']) - max_nt_diff_5p - shift}\t{int(r['end']) + max_nt_diff_3p}\t{r['name']}\t1\t{r['strand']}\n")
            else: 
                bed_file.write(f"{r['chr']}\t{int(r['start']) - max_nt_diff_3p - shift}\t{int(r['end']) + max_nt_diff_5p}\t{r['name']}\t1\t{r['strand']}\n")

def match_chromosomes(path_genomic_file, chr_names):
    # match the exact chr_name in gff 
    for line in fileinput.input(path_genomic_file, inplace = 1): 
        line = line.lower()
        if line.strip().startswith('>'):
            if 'tetraodon8' in line:
                chr = re.search(r">(\w+)", line).group(1)
                print('>' + chr)
            elif chr_names[0] in single_chr_species:
                print('>' + chr_names[0])
            elif 'contig' in chr_names[0]:
                contig_number = int(line.split(' ')[0].split('.')[1])
                contig_name = f"contig{contig_number}"
                if contig_name in chr_names:
                    print('>' + contig_name)
                else:
                    print(line.rstrip())
            else:
                chr = ''
                for chr_name in chr_names:
                    pattern = rf"\b{chr_name}\b"
                    match = re.search(pattern, line)
                    if match: 
                        chr = chr_name
                        break
                if chr == '': 
                    for chr_name in chr_names:
                        if 'chr' in chr_name:
                            number = chr_name.replace('chr', '')
                            pattern = rf"\b(?:chromosome|chr)\s?{number}\b"
                            match = re.search(pattern, line)
                            if match: 
                                chr = chr_name
                                break
                    if chr == '':
                        print(line.rstrip())
                    else:
                        print('>' + chr)
                else:
                    print('>' + chr)
        else: 
            print(line.rstrip())
    
def get_extended_miRNA_coordinates(is_mirbase_gff, match_chr_names, max_nt_diff_5p, max_nt_diff_3p, path_precursors_output_folder, path_genomic_file, path_coords_file):
    chr_names = []

    # read the annotation file from mirbase
    if is_mirbase_gff == True: 
        for line in fileinput.input(path_coords_file, inplace = 1):
            if line.strip().startswith("#"):
                print("")
            else:
                print(line.rstrip())
        genomic_coordinates = pd.read_csv(path_coords_file, sep='\t', header=None, names=["chr", "unknown1", "type", "start", "end", "unknown2", "strand", "unknown3", "details"])
        genomic_coordinates[['id', 'name', 'derives_from']] = genomic_coordinates.apply(lambda r: extract_miR_info(r), axis = 1)
        genomic_coordinates = genomic_coordinates[(genomic_coordinates['type'] == 'miRNA')]
        genomic_coordinates = genomic_coordinates[['chr', 'name', 'start', 'end', 'strand']]
    else:
        genomic_coordinates = pd.read_excel(path_coords_file)

    # Extract chr names 
    genomic_coordinates['chr'] = genomic_coordinates['chr'].str.lower()
    chr_names = list(genomic_coordinates['chr'].unique())
    print(chr_names)

    # create bed coords 
    create_bed_coordinates(genomic_coordinates, max_nt_diff_5p, max_nt_diff_3p, path_precursors_output_folder)   

    # match chr names in genomic coordinate file with genome fasta file 
    if match_chr_names:
        match_chromosomes(path_genomic_file, chr_names)
    # run BEDTools to get extended precursor from genomic coordinates
    abs_path_to_genome_file = os.path.abspath(path_genomic_file)
    abs_path_to_extended_coords_bed_file = os.path.abspath(f'{path_precursors_output_folder}/extended_precursor_coords.bed') 
    abs_path_to_extracted_precursor_seqs_file = os.path.abspath(f'{path_precursors_output_folder}/extended_precursor_seqs.fa')

    # Check if the index file exists
    if os.path.isfile(os.path.abspath(f'{path_genomic_file}.fai')):
    # If it exists, delete the file
        os.remove(os.path.abspath(f'{path_genomic_file}.fai'))

    cmd = f'bedtools getfasta -fi {abs_path_to_genome_file} -bed {abs_path_to_extended_coords_bed_file} -fo {abs_path_to_extracted_precursor_seqs_file} -s -nameOnly'
    subprocess.call(cmd, shell=True)
    # From the output file of BEDTools, save all miRNA names and their extended precursor sequences in a csv file 
    mir_names = [] 
    extended_precursor_seqs = []

    for line in open(abs_path_to_extracted_precursor_seqs_file).readlines():
        if line.strip().startswith('>'):
            mir_names.append(re.search(r'>(.*?)\(', line).group(1))
        else:
            extended_precursor_seq = line.upper()
            extended_precursor_seq = extended_precursor_seq.replace('T', 'U').strip()
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
# Is user using mirbase gff file 
is_mirbase_gff = True if sys.argv[5] == 'True' else False
# Whether user wants to match chromesome of genome file and gff or not
match_chr_names = True if sys.argv[6] == 'True' else False

# Max nt difference at 5p and 3p 
max_nt_diff_5p, max_nt_diff_3p = 0, 0
# create folder if not exists
if not os.path.exists(path_precursors_output_folder):
    os.makedirs(path_precursors_output_folder)

# # Loop through each group
# group_folders = os.listdir(path_summarised_output_folder)
# for group in group_folders:
#     # Get the list of replicate files
#     rep_files = os.listdir(f'{path_summarised_output_folder}/{group}')
#     # Loop through each replicate of that group 
#     for rep_file in rep_files:
#         # Get the replicate name 
#         rep_name = rep_file.split('.')[0]
#         # Read summarised isomiRs file of that replicate 
#         summarised_isomiRs = pd.read_csv(f'{path_summarised_output_folder}/{group}/{rep_file}', encoding='latin-1')
#         # Update max nt difference at 5p and 3p if necessary
#         if max(summarised_isomiRs['5p_nt_diff']) > max_nt_diff_5p: 
#             max_nt_diff_5p = max(summarised_isomiRs['5p_nt_diff'])
#         if max(summarised_isomiRs['3p_nt_diff']) > max_nt_diff_3p:
#             max_nt_diff_3p = max(summarised_isomiRs['3p_nt_diff'])
get_extended_miRNA_coordinates(is_mirbase_gff, match_chr_names, max_nt_diff_5p, max_nt_diff_3p, path_precursors_output_folder, path_genomic_file, path_coords_file)

