#combines results of GWAS from differnet chromosomes into a single input file for Manhattan & QQ plots
#filters out bad snps
import argparse
from os import listdir
from os.path import isfile, join

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--exclude_list")
    parser.add_argument("--result_dir")
    parser.add_argument("--tokeep",nargs="+")
    parser.add_argument("--ending") 
    parser.add_argument("--outf")
    return parser.parse_args()

def get_exclude_dict(fname):
    data=open(fname,'r').read().strip().split('\n')
    exclude_dict=dict()
    for line in data:
        exclude_dict[line]=1
    del data 
    return exclude_dict

def main():
    args=parse_args()
    outf=open(args.outf,'w')
    outf.write('\t'.join(args.tokeep)+'\n')
    print("building exclude dict") 
    exclude_dict=get_exclude_dict(args.exclude_list)
    print('done') 
    files=[f for f in listdir(args.result_dir) if isfile(join(args.result_dir,f))]
    for f in files:
        if f.endswith(args.ending):
            print(str(f))
            data=open(args.result_dir+'/'+f).read().strip().split('\n')
            header=data[0].split()
            keep_indices=[header.index(i) for i in args.tokeep]
            snp_index=header.index('SNP') 
            for line in data[1::]:
                tokens=line.split()
                try:
                    snp=tokens[snp_index]
                    if snp in exclude_dict:
                        continue
                    outf.write('\t'.join([tokens[i] for i in keep_indices])+'\n')
                except:
                    print("failed!:"+str(tokens))
            print('done') 
    
if __name__=="__main__":
    main()
