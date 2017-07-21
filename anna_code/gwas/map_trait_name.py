import sys 
traits=open('ukb_field.tsv','r').read().strip().split('\n')
trait_name_dict=dict()
for line in traits:
    tokens=line.split('\t')
    trait_name_dict[tokens[1]]=tokens[3]
print("made trait dict")
snp_dict=dict()
trait_dict=dict()

data=open(sys.argv[1],'r').read().strip().split('\n')
outf1=open(sys.argv[1]+".named",'w')
outf2=open(sys.argv[1]+".group.snp",'w')
outf3=open(sys.argv[1]+".group.trait","w")
outf1.write(data[0]+'\n')
for line in data[1::]:
    tokens=line.split('\t')
    trait_id=tokens[0]
    trait_name=trait_name_dict[trait_id]
    outf1.write(trait_name+'\t'+'\t'.join(tokens[1::])+'\n')
    chrom=tokens[1]
    snp_name=tokens[2] 
    bp=tokens[3]
    beta=tokens[5]
    pval=tokens[9]
    if trait_name not in trait_dict:
        trait_dict[trait_name]=[]
    trait_dict[trait_name].append(chrom+","+bp+","+beta+","+pval)
    pos_id=chrom+"_"+bp
    if pos_id not in snp_dict:
        snp_dict[pos_id]=[]
    snp_dict[pos_id].append(trait_name+","+beta+","+pval)
print("parsed data")
outf2.write("SNP"+'\t'+"(trait,beta,BonferroniPvalue)"+'\n')
outf3.write("Trait"+'\t'+"(chrom,base,beta,BonferroniPvalue)"+'\n')
for pos_id in snp_dict:
    outf2.write(pos_id+'\t'+'\t'.join(snp_dict[pos_id])+'\n')
for trait_name in trait_dict:
    outf3.write(trait_name+'\t'+'\t'.join(trait_dict[trait_name])+'\n')

