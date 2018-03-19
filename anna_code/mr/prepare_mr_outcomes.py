death_outcome=open('death_status.txt','r').read().strip().split('\n')
cvd_outcome=open('icd_phenotypes_binary.txt','r').read().strip().split('\n')
age_assessed=open('age_attended_assessment_center.txt','r').read().strip().split('\n')

#Aggregate MR outcomes into a single dataframe. 
outf=open('mr_outcome.txt','w')
outf.write('Subject\tDeathStatus\tDeathAge\tCvdStatus\n')
outcomes=dict() 
for line in death_outcome[1::]:
    tokens=line.split('\t')
    fid=tokens[0]
    death_outcomes=set([float(i) for i in tokens[1::]])
    if max(death_outcomes)==-1:
        death_status=0
        death_age=-1 
    else:
        death_status=1
        death_age=max(death_outcomes)
    if fid not in outcomes:
        outcomes[fid]=dict()
    outcomes[fid]['death_status']=str(death_status)
    outcomes[fid]['death_age']=str(death_age)
for line in cvd_outcome[1::]:
    tokens=line.split('\t')
    fid=tokens[0]
    status=int(tokens[-1])-1
    if fid not in outcomes:
        outcomes[fid]=dict()
    outcomes[fid]['cvd_status']=str(status)
for line in age_assessed[1::]:
    tokens=line.split('\t')
    fid=tokens[0]
    age=list(set([float(i) for i in tokens[1::] if i!="NA"]))
    if (len(age)==0):
        age="NA"
    else:
        age=str(max(age))
    if fid not in outcomes:
        outcomes[fid]=dict()
    outcomes[fid]['age_assessed']=str(age) 
    
for fid in outcomes:
    outf.write(fid)
    if 'death_status' in outcomes[fid]:
        outf.write('\t'+outcomes[fid]['death_status'])
    else:
        outf.write('\t')
    if 'death_age' in outcomes[fid]:
        outf.write('\t'+outcomes[fid]['death_age'])
    else:
        outf.write('\t')
    if 'cvd_status' in outcomes[fid]:
        outf.write('\t'+outcomes[fid]['cvd_status'])
    else:
        outf.write('\t')
    if 'age_assessed' in outcomes[fid]:
        outf.write('\t'+outcomes[fid]['age_assessed'])
    else:
        outf.write('\t') 
    outf.write('\n')

        

