#combines results of GWAS from differnet chromosomes into a single input file for Manhattan & QQ plots
#filters out bad snps
import argparse
import pickle
from os import listdir
from os.path import isfile, join

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--exclude_pickle")
    parser.add_argument("--result_dir")
    parser.add_argument("--tokeep",nargs="+")
    parser.add_argument("--ending") 
    parser.add_argument("--outf")
    parser.add_argument("--maf_filter",type=float,default=None)
    parser.add_argument("--maf_pickle_prefix",default=None) 
    return parser.parse_args()

#load from a pickle 
def get_exclude_dict(fname):
    return pickle.load(open(fname,'rb'))

#load from pickles
def get_maf(prefix):
    maf_dict=dict() 
    for chrom in range(1,23):
        chrom=str(chrom)
        maf_dict[chrom]=pickle.load(open(prefix+chrom+".p",'rb'))
    return maf_dict

def main():
    args=parse_args()
    outf=open(args.outf,'w')
    outf.write('\t'.join(args.tokeep)+'\n')
    print("building exclude dict") 
    exclude_dict=get_exclude_dict(args.exclude_pickle)
    print('done')
    if args.maf_pickle_prefix!=None:
        print("building maf dict")
        maf_dict=get_maf(args.maf_pickle_prefix)
        print('done') 
    files=[f for f in listdir(args.result_dir) if isfile(join(args.result_dir,f))]
    if args.maf_filter!=None:
        maf_filter=float(args.maf_filter)    
    for f in files:
        if f.endswith(args.ending):
            print(str(f))
            data=open(args.result_dir+'/'+f).read().strip().split('\n')
            header=data[0].split()
            keep_indices=[header.index(i) for i in args.tokeep]
            snp_index=header.index('SNP')
            chrom_index=header.index('CHR')            
            for line in data[1::]:
                tokens=line.split()
                try:
                    write=True
                    snp=tokens[snp_index]
                    if snp in exclude_dict:
                        write=False
                    if args.maf_filter!=None:
                        chrom=tokens[chrom_index] 
                        maf=float(maf_dict[chrom][snp])
                        if maf<maf_filter:
                            write=False
                    if write==True:
                        outf.write('\t'.join([tokens[i] for i in keep_indices])+'\n')
                except:
                    print("failed!:"+str(tokens))
            print('done') 
    
if __name__=="__main__":
    main()
