#collapses main cause of death to the major letter code for use in Cox Regression Analysis 
data=open("primary_cause_of_death.txt",'r').read().strip().split('\n') 
outf=open("primary_cause_of_death_collapsed.txt",'w')
outf.write(data[0]+'\n')
for line in data[1::]: 
    tokens=line.split('\t') 
    if tokens[-1]!="NA": 
        tokens[-1]=tokens[-1][0] 
    outf.write('\t'.join(tokens)+'\n')

