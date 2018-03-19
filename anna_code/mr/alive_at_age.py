source_data=open("mr_outcome.txt",'r').read().strip().split('\n')
cutoff_age=70
cutoff_age=65
outf=open('mr_outcomes_extra.txt','w')
outf.write(source_data[0]+'\t'+'AliveAt70'+'\t'+'AliveAt65'+'\n')
for line in source_data[1::]:
    tokens=line.split('\t')
    fid=tokens[0]
    
    death_age=None
    known_age=None 
    #is the subject dead or alive?
    if tokens[1]=="1":
        death_age=float(tokens[2])
    else:
        known_age=float(tokens[-1])
    #is the subject older than 70?
    aliveAt70="NA"
    aliveAt65="NA"
    if death_age!=None:
        if death_age <65:
            aliveAt65=0
        if death_age <70:
            aliveAt70=0
    if known_age!=None:
        if known_age>=65:
            aliveAt65=1
        if known_age >=70:
            aliveAt70=1
    outf.write(line+'\t'+str(aliveAt70)+'\t'+str(aliveAt65)+'\n')
    
