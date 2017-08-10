import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="format GWAS hits for sharing w/ 23&Me")
    parser.add_argument('--features_to_keep')
    parser.add_argument('--snp_hits')
    parser.add_argument('--outf')
    return parser.parse_args()

def main():
    args=parse_args()
    hits=open(args.snp_hits,'r').read().strip().split('\n')
    features_to_keep=open(args.features_to_keep,'r').read().strip().split('\n')
    feat_dict=dict()
    for feature in features_to_keep:
        feat_dict[feature]=1
    header=hits[0].split('\t')
    chrom_index=0
    snp_index=1
    pos_index=2
    a1_index=3
    a2_index=4
    maf_index=8
    effect_indices=[] 
    for i in range(8,len(header)):
        cur_label=header[i]
        if cur_label in feat_dict:
            print(str(cur_label))+'-->'+str(i) 
            effect_indices.append(i)
    
    outf=open(args.outf,'w')
    outf.write('SNP\tCHROM\tPOS\tA1\tA2\tMAF\tPHENOTYPE\tP-VALUE\tEFFECT\tSTAT\n')
    print(str(effect_indices))
    for line in hits[1::]:
        tokens=line.split('\t')
        #print(str(tokens))
        chrom=tokens[chrom_index]
        snp=tokens[snp_index]
        pos=tokens[pos_index]
        a1=tokens[a1_index]
        a2=tokens[a2_index]
        maf=tokens[maf_index]
        for effect_index in effect_indices:
            try:
                if tokens[effect_index]!="":
                    effect=tokens[effect_index]
                    stat=tokens[effect_index+1]
                    pval=tokens[effect_index+2]
                    phenotype=header[effect_index]
                    outf.write(snp+'\t'+chrom+'\t'+pos+'\t'+a1+'\t'+a2+'\t'+maf+'\t'+phenotype+'\t'+pval+'\t'+effect+'\t'+stat+'\n')
            except:
                print(str(tokens))
        

if __name__=="__main__":
    main()
    
