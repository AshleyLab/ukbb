original=open("/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/ukb2228.fam",'r').read().strip().split('\n')
samples=open("/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/samples",'r').read().strip().split('\n')
outf=open("/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/ukb2228.filtered.fam",'w')
sample_dict=dict()
for sample in samples:
    sample_dict[sample]=1
for line in original:
    tokens=line.split()
    if tokens[0] in sample_dict:
        outf.write(line+'\n')
        

