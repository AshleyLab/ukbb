#add info value for the key snp's
#optional filter on INFO param
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="add info value to hit snp's, potentially filter on INFO thresh")
    parser.add_argument('--top_snps')
    parser.add_argument('--mfi_prefix')
    parser.add_argument('--outf')
    parser.add_argument('--thresh',type=float,default=0)
    return parser.parse_args()

def main():
    args=parse_args()
    top_snp_data=open(args.top_snps,'r').read().strip().split('\n')
    top_snp_dict=dict()
    for line in top_snp_data[1::]:
        tokens=line.split('\t')
        chrom='chr'+tokens[0]
        pos=tokens[2]
        entry=tuple([chrom,pos])
        top_snp_dict[entry]=line
    print("made top snp dict") 
    outf=open(args.outf,'w')
    outf.write(top_snp_data[0]+'\t'+'Info'+'\n')
    for chrom in range(1,23):
        print(str(chrom))
        chrom='chr'+str(chrom) 
        mfi=open(args.mfi_prefix+str(chrom)+'_v2.txt','r')
        cur_line=mfi.readline()
        while cur_line!="":
            tokens=cur_line.split('\t')
            info=tokens[-1]
            pos=tokens[1]
            entry=tuple([chrom,pos])
            if entry in top_snp_dict:                
                cur_data=top_snp_dict[entry]
                outf.write(cur_data+'\t'+info+'\n')
            cur_line=mfi.readline()
            
    
    

if __name__=="__main__":
    main()
    
