batch_list=open('batches','r').read().strip().split('\n')
batch_dict=dict()
for i in range(len(batch_list)):
    batch_dict[batch_list[i]]=str(i)
original=open('/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/ukb2228.tab.fam','r').read().strip().split('\n')
outf=open('batches.recoded','w')
outf.write(original[0]+'\n')
for line in original[1::]:
    tokens=line.split('\t')
    batch=tokens[-1]
    batch_code=batch_dict[batch]
    tokens[-1]=batch_code
    outf.write('\t'.join(tokens)+'\n')
    
