#combines results of GWAS from differnet chromosomes into a single input file for Manhattan & QQ plots
#filters out bad snps
import argparse
import pickle
from os import listdir
from os.path import isfile, join

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--maf_filter",type=float,default=None)
    parser.add_argument("--maf_pickle_prefix",default=None)
    parser.add_argument("--inputf")
    parser.add_argument("--outf") 
    return parser.parse_args()


#load from pickles
def get_maf(prefix):
    maf_dict=dict() 
    for chrom in range(1,23):
        chrom=str(chrom)
        print(chrom)
        maf_dict[chrom]=pickle.load(open(prefix+chrom+".p",'rb'))
    return maf_dict

def main():
    args=parse_args()
    outf=open(args.outf,'w')
    if args.maf_pickle_prefix!=None:
        print("building maf dict")
        maf_dict=get_maf(args.maf_pickle_prefix)
        print('done')
    data=open(args.inputf,'r').read().strip().split('\n')
    outf.write(data[0]+'\n')
    maf_thresh=float(args.maf_filter) 
    for line in data[1::]:
        tokens=line.split('\t')
        snp=tokens[0]
        chrom=tokens[1]
        maf=float(maf_dict[chrom][snp])
        if maf < maf_thresh:
            continue
        else:
            outf.write(line+'\n')
            
    
if __name__=="__main__":
    main()
