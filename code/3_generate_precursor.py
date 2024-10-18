import pandas as pd 
import subprocess
import re
import os
import re
import sys
import fileinput

def extract_miRBase_miR_info(row):
    """Extract 'id', 'name', 'derives_from' of an miRNA from mirBase gff3. 

    Parameters
    ----------
    row : pd.Series 
        A 1D array that corresponds to a row in mirBase gff3.
    
    Returns 
    -------
    pd.Series 
        Values of 'id', 'name', 'derives_from' 
    """
    details = row['details'].split(';')
    if (row['type'] == "miRNA_primary_transcript"):
        return pd.Series([details[1].split('=')[1], details[2].split('=')[1], ''])
    else:
        return pd.Series([details[1].split('=')[1], details[2].split('=')[1], details[3].split('=')[1]])
    
def create_bed_coordinates(genomic_coordinates, max_nt_diff_5p, max_nt_diff_3p, path_precursor_output_folder):
    """Create a bed file that stores coordinates of extended precursors. 
    
    Parameters
    ----------
    genomic_coordinates : pd.DataFrame 
        A dataframe that stores coordinates of miRNAs. 
    max_nt_diff_5p : int 
        The maximum number of nucleotide difference at 5' end across all isomiRs.
    max_nt_diff_3p : int 
        The maximum number of nucleotide difference at 3' end across all isomiRs.
    path_precursor_output_folder : str 
        Path to the output file.

    Returns 
    -------
    None. A bed file is generated and saved to path_precursor_output_folder.
    
    """
    shift = 1
    with open(f'{path_precursor_output_folder}/extended_precursor_coords.bed', "w+") as bed_file: 
        for _, r in genomic_coordinates.iterrows(): 
            if r['strand'] == '+':
                bed_file.write(f"{r['chr']}\t{int(r['start']) - max_nt_diff_5p - shift}\t{int(r['end']) + max_nt_diff_3p}\t{r['name']}\t1\t{r['strand']}\n")
            else: 
                bed_file.write(f"{r['chr']}\t{int(r['start']) - max_nt_diff_3p - shift}\t{int(r['end']) + max_nt_diff_5p}\t{r['name']}\t1\t{r['strand']}\n")

def match_chromosomes(path_genomic_file, chr_names):
    """Match chromosome names in miRBase or customed coordinates file with those in a given genome file. 

    Parameters
    ----------
    path_genomic_file : str 
        Path to the genome file. 
    chr_names : list 
        List of chromosome names in miRBase or customed coordinates file.

    Returns
    -------
    None. The genome file's chromosome names are replaced with those from miRBase or customed coordinates file (if they match).
    """
    for line in fileinput.input(path_genomic_file, inplace = 1): 
        line = line.lower()
        if line.strip().startswith('>'):
            if 'tetraodon8' in line:
                chr = re.search(r">(\w+)", line).group(1)
                print('>' + chr)
            elif len(chr_names) == 1:
                print('>' + chr_names[0])
            elif 'contig' in chr_names[0]:
                contig_number = int(re.search(r'\.([0-9]+)\s', line).group(1))
                contig_name = f"contig{contig_number}"
                if contig_name in chr_names:
                    print('>' + contig_name)
                else:
                    print(line.rstrip())
            else:
                chr = ''
                for chr_name in chr_names:
                    match = re.search(rf"\b{chr_name}\b", line)
                    if match: 
                        chr = chr_name
                        break
                if chr == '': 
                    for chr_name in chr_names:
                        if 'chr' in chr_name:
                            number = chr_name.replace('chr', '')
                            match = re.search(rf"\b(?:chromosome|chr)\s?{number}\b", line)
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
    """Extract the extended precursor sequences for miRNAs.

    Parameters
    ----------
    is_mirbase_gff : boolean
        Is the miRNA coordinates from miRBase gff file ? 
    match_chr_names : boolean 
        Is matching chromosome names required ? It must be true if the chromosome names in genome file are not identical to those in miR coordinates file. 
    max_nt_diff_5p : int
        The maximum number of nucleotide difference at 5' end across all isomiRs.
    max_nt_diff_3p : int
        The maximum number of nucleotide difference at 3' end across all isomiRs.
    path_precursors_output_folder : str 
        Path to the folder that stores extracted precursor sequences outputs. 
    path_genomic_file : str 
        Path to the genome file. 
    path_coords_file : str 
        Path to the miRNA coordinates file. 

    Returns 
    -------
    None. All extracted sequences are saved in a csv file in the path_precursors_output_folder.

    """
    chr_names = []
    genomic_coordinates = pd.DataFrame()

    if is_mirbase_gff == True: 
        for line in fileinput.input(path_coords_file, inplace = 1):
            if line.strip().startswith("#"):
                print("")
            else:
                print(line.rstrip())
        genomic_coordinates = pd.read_csv(path_coords_file, sep='\t', header=None, names=["chr", "unknown1", "type", "start", "end", "unknown2", "strand", "unknown3", "details"])
        genomic_coordinates[['id', 'name', 'derives_from']] = genomic_coordinates.apply(lambda r: extract_miRBase_miR_info(r), axis = 1)
        genomic_coordinates = genomic_coordinates[(genomic_coordinates['type'] == 'miRNA')]
        genomic_coordinates = genomic_coordinates[['chr', 'name', 'start', 'end', 'strand']]
    else:
        genomic_coordinates = pd.read_excel(path_coords_file)

    # Match chr names in genomic coordinate file with those in genome fasta file 
    if match_chr_names:
        genomic_coordinates['chr'] = genomic_coordinates['chr'].astype(str)
        genomic_coordinates['chr'] = genomic_coordinates['chr'].str.lower()        
        chr_names = list(genomic_coordinates['chr'].unique())
        match_chromosomes(path_genomic_file, chr_names)
    
    # Create a bed file that stores coordinates of extended precursors 
    create_bed_coordinates(genomic_coordinates, max_nt_diff_5p, max_nt_diff_3p, path_precursors_output_folder)   

    abs_path_to_genome_file = os.path.abspath(path_genomic_file)
    abs_path_to_extended_coords_bed_file = os.path.abspath(f'{path_precursors_output_folder}/extended_precursor_coords.bed') 
    abs_path_to_extracted_precursor_seqs_file = os.path.abspath(f'{path_precursors_output_folder}/extended_precursor_seqs.fa')

    # Check if the index file exists
    if os.path.isfile(os.path.abspath(f'{path_genomic_file}.fai')):
    # If it exists, delete the file
        os.remove(os.path.abspath(f'{path_genomic_file}.fai'))

    # Run BEDTools to get extended precursor from genomic coordinates
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
            miRNA_seq = extended_precursor_seq[max_nt_diff_5p:len(extended_precursor_seq) - max_nt_diff_3p]
            extended_precursor_seq = extended_precursor_seq.replace(miRNA_seq, miRNA_seq.lower())
            extended_precursor_seqs.append(extended_precursor_seq)
    pd.DataFrame({'mir_name': mir_names, 'extended_precursor_seq': extended_precursor_seqs}).to_csv(f'{path_precursors_output_folder}/{max_nt_diff_5p}_{max_nt_diff_3p}_extended_precursor_seqs.csv', index=False)

print("Running 3_generate_precursor.py script...")

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
# Does user want to match chromesome of genome file and gff or not
match_chr_names = True if sys.argv[6] == 'True' else False

# Max nt difference at 5p and 3p 
max_nt_diff_5p, max_nt_diff_3p = 0, 0
# Create folder if not exists
if not os.path.exists(path_precursors_output_folder):
    os.makedirs(path_precursors_output_folder)

# List of group folders 
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
get_extended_miRNA_coordinates(is_mirbase_gff, match_chr_names, max_nt_diff_5p, max_nt_diff_3p, path_precursors_output_folder, path_genomic_file, path_coords_file)

