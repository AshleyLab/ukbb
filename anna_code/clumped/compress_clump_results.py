#compresses clumped SNPs for a particular feature set to a single matrix
import argparse
from os import listdir
from os.path import isfile, join
import math

def parse_args():
    parser=argparse.ArgumentParser(description="clump features into a single file")
    parser.add_argument("--source_dir",default="clumped")
    parser.add_argument("--outf",default="clumped.matrix.txt")
    return parser.parse_args() 

def main():
    args=parse_args()
    snp_dict=dict()
    mypath=args.source_dir
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    snp_files=[]
    for f in onlyfiles:
        if f.endswith('clumped'):
            snp_files.append(f)
    outf=open(args.outf,'w')
    outf.write('SNP\t'+'\t'.join(snp_files)+'\n')
    for f in snp_files:
        data=open(args.source_dir+'/'+f,'r').read().strip().split('\n')
        for line in data:
            tokens=line.split('\t')
            #print(str(tokens))
            if len(tokens)<2:
                continue 
            snp=tokens[0]
            pval=tokens[1]
            if snp not in snp_dict:
                snp_dict[snp]=dict()
            snp_dict[snp][f]=pval
    for snp in snp_dict:
        outf.write(snp)
        for f in snp_files:
            if f in snp_dict[snp]:
                outf.write('\t'+str(-10*math.log10(float(snp_dict[snp][f]))))
            else:
                outf.write('\t')
        outf.write('\n')
    
if __name__=="__main__":
    main()
    
