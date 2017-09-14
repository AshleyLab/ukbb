#filter PLINK hits by maf cutoff
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="filter PLINK results by MAF")
    parser.add_argument("--input")
    parser.add_argument("--thresh",type=float,default=5.0e-8)
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    data=open(args.input,'r').read().strip().split('\n')
    outf=open(args.outf,'w')
    header=data[0].split()
    outf.write('\t'.join(header))
    for line in data[1::]:
        tokens=line.split()
        pval=tokens[-1]
        if pval=="NA":
            continue
        pval=float(pval)
        print(str(pval))
        if pval <= args.thresh:
            outf.write('\t'.join(tokens)+'\n')
            
        
if __name__=="__main__":
    main()
    
