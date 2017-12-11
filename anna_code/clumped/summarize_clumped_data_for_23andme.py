import argparse
import os
from os import listdir
from os.path import isfile, join

def parse_args():
    parser=argparse.ArgumentParser(description="annotate clumped GWAS summaries for reporting to 23&Me")
    parser.add_argument("--files_to_summarize")
    parser.add_argument("--maf_base_dir")
    parser.add_argument("--gwas_dirs")
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    outf=open(args.outf,'w') 
    files_to_summarize=open(args.files_to_summarize,'r').read().strip().split('\n')
    gwas_dirs=open(args.gwas_dirs,'r').read().strip().split('\n') 
    snp_dict=dict()
    for i in range(len(files_to_summarize)):
        print(str(i))
        f=files_to_summarize[i] 
        sig_snps=open(f,'r').read().strip().split('\n')
        for line in sig_snps:
            tokens=line.split('\t')
            try:
                rsid=tokens[0]
                pval=tokens[1]
                if rsid not in snp_dict:
                    snp_dict[rsid]=dict()
                snp_dict[rsid][i]=1
            except:
                continue 
    print("generated SNP dict") 
    #get the GWAS statistics for SNP/feature combination
    for i in range(len(gwas_dirs)):
        print(str(i)) 
        cur_dir=gwas_dirs[i] 
        onlyfiles = [f for f in listdir(cur_dir) if isfile(join(cur_dir, f))]
        for f in onlyfiles:
            if ((f.endswith("linear")) or (f.endswith("logistic"))):
                print(str(f))
                data=open(cur_dir+'/'+f,'r').read().strip().split("\n")
                for line in data[1::]:
                    tokens=line.split()
                    rsid=tokens[2]
                    if ((rsid in snp_dict) and (i in snp_dict[rsid])):
                        #record the entry
                        snp_dict[rsid][i]=line
    print("added GWAS hit information for significant SNP hits")
    #get the minor allele frequency information for the significant snps.
    for chrom in range(1,23):
        print(str(chrom)) 
        cur_maf=open(args.maf_base_dir+'/ukb_mfi_chr'+str(chrom)+"_v2.txt",'r')
        for line in cur_maf:
            tokens=line.split()
            if tokens[0] in snp_dict:
                rsid=tokens[0]
                for i in snp_dict[rsid]:
                    outf.write(line+'\t'+files_to_summarize[i]+'\t'+snp_dict[rsid][i]+'\n')
    
            
    
if __name__=="__main__":
    main()
    
