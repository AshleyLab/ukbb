#generates pickles of
# 1) exclude snps
# 2) maf
# 3) GWAS catalog
#import cPickle as pickle
import pickle
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="helper functions to generate useful data pickles")
    parser.add_argument("--maf_base")
    parser.add_argument("--chrom")
    parser.add_argument("--outf")
    return parser.parse_args()

def get_exclude_dict(fname):
    exclude_dict=dict()
    data=open(fname,'r').read().strip().split('\n')
    for line in data:
        exclude_dict[line]=1
    del data
    return exclude_dict

def get_maf_dict(fname):
    maf_dict=dict()
    for chrom in range(1,23):
        print("maf for chrom:"+str(chrom))
        data=open(fname+str(chrom)+"_v2.txt",'r').read().strip().split('\n')
        for line in data:
            tokens=line.split()
            if len(tokens)<1:
                print(str(tokens))
            else:
                snp=tokens[0]
                maf=tokens[-2]
                maf_dict[snp]=maf
    del data
    return maf_dict 


def get_maf_dict_chrom(fname):
    maf_dict=dict()
    data=open(fname,'r').read().strip().split('\n')
    for line in data:
        tokens=line.split()
        if len(tokens)<1:
            continue
        else:
            snp=tokens[0]
            maf=tokens[-2]
            maf_dict[snp]=maf
    del data
    return maf_dict 

def get_gc_dict(fname,exclude_snps):
    data=open(fname,'r').read().strip().split('\n')
    gc_dict=dict()
    for line in data[1::]:
        tokens=line.split('\t')
        snp=tokens[2]
        if snp not in exclude_snps: 
            pmid=tokens[-3]
            trait=tokens[-2]
            gc_dict[snp]=[pmid,trait]
    del data
    return gc_dict

def main():
    args=parse_args()
    maf_base=args.maf_base
    chrom=args.chrom
    fname=maf_base+chrom+"_v2.txt"
    maf_dict=get_maf_dict_chrom(fname)
    pickle.dump(maf_dict,open(args.outf,'wb'))
    
#exclude_dict=get_exclude_dict("/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/hrc/bad.snps.txt")
#pickle.dump(exclude_dict,open("exclude_dict.p","wb"))
#maf_dict=get_maf_dict("/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/EGAD00010001225/001/ukb_mfi_chr")
#print("writing pickle!") 
#pickle.dump(maf_dict,open("maf_dict.p","wb"))
#gwas_dict=get_gc_dict("/scratch/PI/euan/common/gwascatalog/gwas_catalog_may2017.slim.txt",exclude_dict)
#pickle.dump(gwas_dict,open("gwas_dict.p","wb"))


if __name__=="__main__":
    main()
    
