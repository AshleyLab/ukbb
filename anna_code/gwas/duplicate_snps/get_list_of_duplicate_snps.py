#lists snps that are duplicated in the batch 2 bim files
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="lists snps that are duplicated in the batch 2 bim files")
    parser.add_argument("--source_prefix",default="/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/plink/ukb_imp_chr")
    parser.add_argument("--source_suffix",default="_v2.bim") 
    parser.add_argument("--outf",default="duplicate_snps.txt")
    parser.add_argument("--chrom") 
    return parser.parse_args()

def main():
    args=parse_args()
    dups=dict()
    chrom=args.chrom
    outf=open(args.outf,'w') 
    snps=open(args.source_prefix+chrom+args.source_suffix,'r').read().strip().split('\n')
    for line in snps:
        cursnp=line.split()[1]
        if cursnp in dups:
            outf.write(cursnp+'\n')
        else:
            dups[cursnp]=1
    
            
if __name__=="__main__":
    main()
    
