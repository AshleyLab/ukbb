import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="classify each hit SNP as called or imputed")
    parser.add_argument("--hits")
    parser.add_argument("--chrom_column",type=int)
    parser.add_argument("--pos_column",type=int)
    parser.add_argument("--snp_column",type=int)
    parser.add_argument("--calls_dir")
    parser.add_argument("--outf")
    return parser.parse_args()
def main():
    args=parse_args()
    data=open(args.hits,'r').read().strip().split('\n')
    outf=open(args.outf,'w')
    outf.write('SNP\tChrom\tPos\tStatus\n')
    chrom_col=args.chrom_column
    pos_col=args.pos_column
    snp_col=args.snp_column
    calls=dict()
    for chrom in range(1,23):
        print("parsing calls for chrom:"+str(chrom))
        call_data=open(args.calls_dir+"/"+"ukb_snp_chr"+str(chrom)+"_v2.bim",'r').read().strip().split('\n')
        for line in call_data:
            tokens=line.split('\t')
            chrom=tokens[0]
            pos=tokens[3]
            entry=tuple([chrom,pos])
            calls[entry]=1

    for line in data[1::]:
        tokens=line.split('\t')
        chrom=tokens[chrom_col]
        pos=tokens[pos_col]
        entry=tuple([chrom,pos])
        if entry in calls:
            outf.write(tokens[snp_col]+'\t'+chrom+'\t'+pos+'\t'+'Called'+'\n')
        else:
            outf.write(tokens[snp_col]+'\t'+chrom+'\t'+pos+'\t'+'Imputed'+'\n')
        
if __name__=="__main__":
    main()
    
    
