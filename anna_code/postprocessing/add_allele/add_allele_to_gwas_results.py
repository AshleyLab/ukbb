import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="Add allele info to GWAS hits")
    parser.add_argument('--frequencies')
    parser.add_argument('--gwas_output',help="file containing full path to all gwas outputs you would like analyzed")
    parser.add_argument("--chrom",default=None,type=int)
    return parser.parse_args()

def read_freqs(fname,chrom):
        maf_dict=dict()
        if chrom==None:
            for chrom in range(1,23):
                print("maf for chrom:"+str(chrom))
                data=open(fname+str(chrom)+"_v2.txt",'r').read().split('\n')
                for line in data:
                    tokens=line.split()
                    if len(tokens)<1:
                        print(str(tokens))
                    else:
                        snp=tokens[0]
                        maf_dict[snp]=tokens[2:4]
        else:
            print("maf for chrom:"+str(chrom))
            data=open(fname+str(chrom)+"_v2.txt",'r').read().strip().split('\n')
            for line in data:
                tokens=line.split()
                if len(tokens)<1:
                    print(str(tokens))
                else:
                    snp=tokens[0]
                    maf_dict[snp]=tokens[2:4]
        del data 
        return maf_dict
                                                                                                                                                        
def annotate(fname,freq_dict,chrom):
    print("parsing file:"+str(fname)) 
    data=open(fname,'r').read().strip().split('\n')
    if chrom==None:
        outf=open(fname+".withAlleles",'w')
    else:
        outf=open(fname+'.'+str(chrom)+".withAlleles",'w')
    header=data[0].split()
    snp_index=header.index('SNP')
    chrom_index=header.index('CHR') 
    a2_only=False
    if 'A1' in header:
        a2_only=True
        a1_index=header.index('A1')
        outf.write('\t'.join(header)+'\tA2\n') 
    else: 
        outf.write('\t'.join(header)+'\tA1\tA2\n') 
    for line in data[1::]:
        try:
            tokens=line.split('\t')
            cur_chrom=int(tokens[chrom_index])
            if cur_chrom!=chrom:
                continue 
            cur_freqs=freq_dict[tokens[snp_index]]
            if (a2_only==True):
                if cur_freqs[0]==tokens[a1_index]:
                    a2=cur_freqs[1]
                else:
                    a2=cur_freqs[0] 
                outf.write('\t'.join(tokens)+'\t'+a2+'\n')
            else:
                outf.write('\t'.join(tokens)+'\t'+'\t'.join(cur_freqs)+'\n')        
        except:
            print("WARNING! Skipping line:"+line)
def main():
    args=parse_args()
    #read in the frequency information
    freq_dict=read_freqs(args.frequencies,args.chrom)
    file_list=open(args.gwas_output,'r').read().strip().split('\n')
    for fname in file_list:
        annotate(fname,freq_dict,args.chrom)
        
    
    

if __name__=="__main__":
    main()
    
