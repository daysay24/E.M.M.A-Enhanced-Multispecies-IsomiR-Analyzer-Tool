import pandas as pd
import os



def extract_miR_info(row):
    details = row['details'].split(';')
    if (row['type'] == "miRNA_primary_transcript"):
        return pd.Series([details[1].split('=')[1], details[2].split('=')[1], ''])
    else:
        return pd.Series([details[1].split('=')[1], details[2].split('=')[1], details[3].split('=')[1]])


# read extracted output file 
extended_output = pd.read_csv('./output/3_precursors/0_0_extended_precursor_seqs.csv')
# # read mirbase mir sequences 
# mir_names=[]
# mir_seqs = []
# with open('./mature.fa', 'r') as f: 
#     for line in f: 
#         if '>' in line:
#             mir_names.append(line.split(' ')[0].replace('>', ''))
#         else:
#             mir_seqs.append(line.strip().lower())
# mature = pd.DataFrame({'mir_name': mir_names, 'mir_seq': mir_seqs})
# mature.to_csv('./mature.csv', index=False)

# get speices name 
test_files = os.listdir('./test')
species = [file for file in test_files if 'fna' in file][0].split('/')[-1].split('.')[0]

# Read mature csv 
mature = pd.read_csv('./mature.csv')

# mir_name mature 
mmu_mir_name_1 = set(mature[mature['mir_name'].str.startswith(species)]['mir_name'])
# mir_name output 
mmu_mir_name_2 = set(extended_output['mir_name'])
# mir_name gff
genomic_coordinates = pd.read_csv(f'./test/{species}.gff3', sep='\t', header=None, names=["chr", "unknown1", "type", "start", "end", "unknown2", "strand", "unknown3", "details"])
genomic_coordinates[['id', 'name', 'derives_from']] = genomic_coordinates.apply(lambda r: extract_miR_info(r), axis = 1)
genomic_coordinates = genomic_coordinates[(genomic_coordinates['type'] == 'miRNA')]
genomic_coordinates = genomic_coordinates[['chr', 'name', 'start', 'end', 'strand']]
mmu_mir_name_3 = set(genomic_coordinates['name'])
print("In mature.fa, but not in extracted output: ", mmu_mir_name_1 - mmu_mir_name_2)
print("======================")
print("In mature.fa, but not in gff: ", mmu_mir_name_1 - mmu_mir_name_3)
print("======================")
print("In gff, but not in extracted output: ", mmu_mir_name_2 - mmu_mir_name_3)
print("======================")

# join 2 tables by mir name 
merged = extended_output.merge(mature, how='inner', on='mir_name')
print('Extracted sequences are the same as mature.fa sequences: ', merged['extended_precursor_seq'].equals(merged['mir_seq']))
print("======================")
# compare 2 columns 
for _, r in merged.iterrows():
    if r['extended_precursor_seq'] != r['mir_seq']:
        print("Mirname whose sequence is difference between mature.fa and extracted output: ", r['mir_name'])