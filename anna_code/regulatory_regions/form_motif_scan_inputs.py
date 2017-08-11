import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="formulate inputs for motif scanning")
    parser.add_argument("--snp_hits")
    parser.add_argument("--flank",type=int,default=30)
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    hits=open(args.snp_hits,'r').read().strip().split('\n')
    header=hits[0].split('\t')
    chrom_index=header.index('CHR')
    pos_index=header.index('POS')
    snp_index=header.index('SNP')
    outf=open(args.outf,'w')
    for line in hits[1::]:
        tokens=line.split('\t')
        chrom=tokens[chrom_index]
        pos=int(tokens[pos_index])
        snp=tokens[snp_index]
        startpos=pos-args.flank
        endpos=pos+args.flank
        outf.write(chrom+'\t'+str(startpos)+'\t'+str(endpos)+'\t'+snp+'\n')
        

if __name__=="__main__":
    main()
    
