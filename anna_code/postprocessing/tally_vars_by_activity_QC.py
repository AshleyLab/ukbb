data=open("top_snps_activity_merged.annotations",'r').read().strip().split('\n')
hits=dict()
for line in data[1::]:
    tokens=line.split('\t')
    snp=tokens[0]
    features=tokens[4].split(',')
    for f in features:
        if f not in hits:
            hits[f]=1
        else:
            hits[f]+=1
for key in hits:
    print(str(key)+'\t'+str(hits[key]))
    
