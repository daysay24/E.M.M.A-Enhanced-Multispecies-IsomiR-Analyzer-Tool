import pandas as pd
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

# Read mature csv 
mature = pd.read_csv('./mature.csv')
# join 2 tables by mir name 
merged = extended_output.merge(mature, how='inner', on='mir_name')
print(merged['extended_precursor_seq'].equals(merged['mir_seq']))
# compare 2 columns 
for _, r in merged.iterrows():
    if r['extended_precursor_seq'] != r['mir_seq']:
        print(r['mir_name'])