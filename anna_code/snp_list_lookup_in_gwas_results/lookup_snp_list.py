import argparse 
import pandas as pd 
import pdb 
def parse_args(): 
    parser=argparse.ArgumentParser(description="query a list of variants from GWAS results") 
    parser.add_argument("--var_list") 
    parser.add_argument("--gwas_hits",nargs="+")
    parser.add_argument("--thresh",default=0.05,type=float)
    parser.add_argument("--phenotype") 
    parser.add_argument("--outf") 
    parser.add_argument("--index_col",type=int,default=2)
    return parser.parse_args()

def main(): 
    args=parse_args() 
    variants=[i.split('\t')[0] for i in open(args.var_list,'r').read().strip().split('\n')]
    outf=open(args.outf,'w')
    for filename in args.gwas_hits: 
        data=pd.read_table(filename,index_col=args.index_col,delim_whitespace=True)
        print("checking "+str(filename))
        for variant in variants:
            try:
                pval=data['P'].loc[variant]
                outf.write(variant+'\t'+str(pval)+'\t'+args.phenotype+'\n')
            except: 
                print("skipping:"+str(variant))
if __name__=="__main__": 
    main() 
