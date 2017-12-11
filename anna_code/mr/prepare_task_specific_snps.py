#subsets the mr_snps.txt dataframe to select the snps that are significant for each exposure.
import pandas as pd
import argparse
from os import listdir
from os.path import isfile,join

def parse_args():
    parser=argparse.ArgumentParser("subsets the mr_snps.txt dataframe to select the snps that are significant for each exposure.")
    parser.add_argument("--feature_file_dir")
    parser.add_argument("--snp_mat")
    parser.add_argument("--output_dir") 
    return parser.parse_args()

def main():
    args=parse_args()
    data=pd.read_table(args.snp_mat,sep='\t')
    print('read in data!')
    feature_files=[f for f in listdir(args.feature_file_dir) if isfile(join(args.feature_file_dir,f))]
    for f in feature_files:
        full_path=args.feature_file_dir+'/'+f
        snps=[line.split('\t')[0] for line in open(full_path,'r').read().strip().split('\n')]
        if ((snps[0]=="") and (len(snps)==1)):
            print("skipping:"+str(f))
            continue 
        subset=data[snps]
        outf=args.output_dir+'/'+f.split('.')[0]
        subset.to_csv(outf,sep='\t')
        print('wrote:'+str(outf))
        

if __name__=="__main__":
    main()
    
