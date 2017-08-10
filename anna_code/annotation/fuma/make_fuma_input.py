#generates files in format for FUMA analysis
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="generate input for fuma analysis")
    parser.add_argument("--inputf")
    parser.add_argument("--outf")
    parser.add_argument('--summary_with_a1_a2')
    return parser.parse_args()

def main():
    args=parse_args()
    summary=open(args.summary_with_a1_a2,'r').read().strip().split('\n')
    summary_dict=dict()
    for line in summary[1::]:
        tokens=line.split('\t')
        snp=tokens[1]
        a1=tokens[3]
        a2=tokens[4]
        summary_dict[snp]=[a1,a2]
    print(str(summary_dict))
    data=open(args.inputf,'r').read().strip().split('\n')
    outf=open(args.outf,'w')
    header=data[0].split('\t')
    outf.write('P-value\trsID\tchromosome\tposition\tN\tA1\tA2\n')
    p_index=header.index('pvalue')
    snp_index=header.index('snp')
    chrom_index=header.index('chr')
    pos_index=header.index('pos')
    nmiss_index=header.index('nmiss')
    for line in data[1::]:
        tokens=line.split('\t')
        try:
            a1=summary_dict[tokens[snp_index]][0] 
            a2=summary_dict[tokens[snp_index]][1] 
            outf.write(tokens[p_index]+'\t'+tokens[snp_index]+'\t'+tokens[chrom_index]+'\t'+tokens[pos_index]+'\t'+tokens[nmiss_index]+'\t'+a1+'\t'+a2+'\n')
        except:
            #this will catch all snps in the bad snp exclusion list that are present in the individual top snp files, but absent from the merged file. 
            continue
    
if __name__=="__main__":
    main()
    
