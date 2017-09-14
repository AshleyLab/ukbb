covars=open("covariates.emmi.txt",'r').read().strip().split('\n')
phenotypes=open("ukb_PA_GS_SLEEP.pheno.recoded",'r').read().strip().split('\n')
outf=open("covariates.emmi.filtered.txt",'w')
phenotype_dict=dict()
for line in phenotypes[1::]:
    tokens=line.split('\t')
    iid=tokens[0]
    phenotype_dict[iid]=1
covar_header=covars[0]
outf.write(covar_header+'\n')
for line in covars[1::]:
    tokens=line.split('\t')
    cur_iid=tokens[0]
    if cur_iid in phenotype_dict:
        outf.write(line+'\n')
        
