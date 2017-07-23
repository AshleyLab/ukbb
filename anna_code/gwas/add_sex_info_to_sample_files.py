sample_file=open("/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/ukb2228.sample",'r').read().strip().split('\n') 
fam_file=open("/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/ukb2228.fam",'r').read().strip().split('\n') 
outf=open("sample.augmented.txt",'w')
fam_dict=dict()
for line in fam_file:
    tokens=line.split()
    fam_dict[tokens[0]]=tokens[4]
outf.write(sample_file[0]+'\tsex\n')
for line in sample_file[1::]:
    tokens=line.split()
    subject=tokens[0]
    if subject not in fam_dict:
        sex="0"
    else:
        sex=fam_dict[subject]
    outf.write(line+'\t'+sex+'\n')
    
