import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="Add allele info to GWAS hits")
    parser.add_argument('--frequencies')
    parser.add_argument('--gwas_output',help="file containing full path to all gwas outputs you would like analyzed")
    parser.add_argument('--t',default=1)
    return parser.parse_args()

def read_freqs(freq_file):
    freqs=open(freq_file,'r').read().strip().split("\n")
    print("loaded frequency file") 
    freq_dict=dict()
    for line in freqs:
        tokens=line.split()
        freq_dict[tokens[0]]=tokens[2:4]
    print("constructed frequency dictionary")
    del freqs 
    return freq_dict

def annotate(fname,freq_dict):
    print("parsing file:"+str(fname)) 
    data=open(fname,'r').read().strip().split('\n') 
    outf=open(fname+".withAlleles",'w')
    header=data[0].split()
    snp_index=header.index('SNP')
    a2_only=False 
    if 'A1' in header:
        a2_only=True
        a1_index=header.index('A1')
        outf.write(data[0]+'\tA2\n') 
    else: 
        outf.write(data[0]+'\tA1\tA2\n') 
    for line in data[1::]:
        tokens=line.split()
        cur_freqs=freq_dict[tokens[snp_index]]
        if a2_only:
            a1=tokens[a1_index]
            cur_freqs.remove(a1) 
            outf.write(line+'\t'+'\t'.join(cur_freqs)+'\n')
        else:
            outf.write(line+'\t'+'\t'.join(cur_freqs))        
                   
def main():
    args=parse_args()
    #read in the frequency information
    freq_dict=read_freqs(args.frequencies)
    file_list=open(args.gwas_output,'r').read().strip().split('\n')
    for fname in file_list:
        annotate(fname,freq_dict)
        
    
    

if __name__=="__main__":
    main()
    
