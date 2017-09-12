#converts subject id's in phenotype file from application 13721 to 2228
bridge=open("bridging_file.csv",'r').read().strip().split('\n')
bridge_dict=dict()
for line in bridge[1::]:
    tokens=line.split(',')
    bridge_dict[tokens[1]]=tokens[0]
data=open("/oak/stanford/groups/euan/projects/ukbb/data/phenotype/ukb_PA_GS_SLEEP.pheno",'r').read().strip().split('\n')
outf=open("/oak/stanford/groups/euan/projects/ukbb/data/phenotype/ukb_PA_GS_SLEEP.pheno.recoded",'w')
outf.write(data[0]+'\n')
for line in data[1::]:
    tokens=line.split('\t')
    subject=tokens[0]
    if subject not in bridge_dict:
        print(str(tokens))
        continue 
    recoded_subject=bridge_dict[subject]
    tokens[0]=recoded_subject
    tokens[1]=recoded_subject
    outf.write(' '.join(tokens)+'\n')
    
