#binarizes categorical phenotypes at a provided threshold
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="binarize categorical phenotype at a provided threshold")
    parser.add_argument("--inputf")
    parser.add_argument("--thresh",type=int)
    parser.add_argument("--outf")
    parser.add_argument("--missing",type=int,default=-1000)
    parser.add_argument("--first_data_column",type=int,default=2) 
    return parser.parse_args()

def main():
    args=parse_args()
    data=open(args.inputf,'r').read().strip().split('\n')
    header=data[0]
    outf=open(args.outf,'w') 
    outf.write(header+'\n')
    for line in data[1::]:
        tokens=line.split('\t')
        outstring=tokens[0:args.first_data_column] 
        values=[float(i) for i in tokens[args.first_data_column::]]
        for v in values:
            if v==args.missing:
                outstring.append(str(int(v)))
            elif v > args.thresh:
                outstring.append('1')
            else:
                outstring.append('0')
        outf.write('\t'.join(outstring)+'\n')

if __name__=="__main__":
    main()
    
    
