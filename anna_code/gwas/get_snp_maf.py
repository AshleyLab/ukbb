import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="get population maf values for PLINK qassoc SNP hits")
    parser.add_argument("--plink_result")
    parser.add_argument("--freq_dir")
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    plink_result=open(args.plink_result,'r').read().strip().split('\n')
    snp_dict=dict()
    header=plink_result[0].split('\t')
    print(str(header))
    chrom_index=header.index('CHR')
    snp_index=header.index('SNP')
    for line in plink_result[1::]:
        tokens=line.split('\t')
        chrom=tokens[chrom_index]
        snp=tokens[snp_index]
        if chrom not in snp_dict:
            snp_dict[chrom]=dict()
        snp_dict[chrom][snp]=line
    print ("populated SNP dict from PLINK results")
    outf=open(args.outf,'w')
    outf.write('\t'.join(header)+'\tmaf\n')
    for chrom in snp_dict:
        print("getting maf for chrom:"+str(chrom))
        maf_file=open(args.freq_dir+'/chrom'+chrom+".frq",'r').read().strip().split('\n')
        for line in maf_file[1::]:
            tokens=line.split()
            print(str(tokens))
            snp=tokens[1]
            if snp in snp_dict[chrom]:
                maf=tokens[4] 
                outf.write(snp_dict[chrom][snp]+'\t'+maf+'\n')
        
if __name__=="__main__":
    main() 
