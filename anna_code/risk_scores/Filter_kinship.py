import pandas as pd
import pdb 
kinship = pd.read_csv("/oak/stanford/groups/euan/projects/ukbb/code/anna_code/risk_scores/plink.genome", delim_whitespace=True)
first_degree=kinship[kinship['PI_HAT']>.4]
pdb.set_trace() 
first_degree.to_csv("first_degree_relatives.txt") 

kin_dict={} 
throw_away={} 

for index,row in first_degree.iterrows(): 
    p1=row['IID1']
    p2=row['IID2'] 
    if p1 not in kin_dict: 
        kin_dict[p1]=[p2] 
    else: 
        kin_dict[p1].append(p2) 
    if p2 not in kin_dict: 
        kin_dict[p2]=[p1] 
    else: 
        kin_dict[p2].append(p1) 
print("made kinship dictionary") 

for person in kin_dict: 
    if person not in throw_away: 
        for relative in kin_dict[person]: 
            throw_away[relative]=1 
outf=open("throw_away.tsv",'w') 
throw_away_list=list(throw_away.keys())
outf.write('\n'.join([str(i) for i in throw_away_list]))
