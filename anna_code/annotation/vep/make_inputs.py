#generates inputs for vep analysis
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="generates inputs for vep analysis")
    parser.add_argument("--inputf")
    parser.add_argument("--outputf")
    return parser.parse_args()

def main():
    args=parse_args()
    data=open(args.inputf,'r').read().strip().split('\n')
    header=data[0].split('\t')
    chrom_index=header.index('CHR')
    pos_index=header.index('POS')
    a1_index=header.index('A1')
    a2_index=header.index('A2')
    outf=open(args.outputf,'w')
    for line in data[1::]:
        tokens=line.split('\t')
        print(str(tokens))
        outf.write(tokens[chrom_index]+'\t'+tokens[pos_index]+'\t'+tokens[pos_index]+'\t'+tokens[a1_index]+'/'+tokens[a2_index]+'\t'+'1'+'\n')
if __name__=="__main__":
    main()
