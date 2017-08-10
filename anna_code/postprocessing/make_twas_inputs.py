import argparse
import pdb 
def parse_args():
    parser=argparse.ArgumentParser(description="Add allele info to GWAS hits")
    parser.add_argument('--frequencies')
    parser.add_argument('--gwas_output',help="file containing full path to all gwas outputs you would like analyzed")
    parser.add_argument('--t',default=1)
    parser.add_argument('--keep_columns',nargs="*")
    parser.add_argument('--exclude_list',default=None,help="list of SNPs to exclude from the final report")
    parser.add_argument('--outdir',help="output directory to store the twas inputs")
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

def annotate(fname,freq_dict,keep_columns,exclude_dict,outdir):
    print("parsing file:"+str(fname)) 
    data=open(fname,'r').read().strip().split('\n') 
    outf=open(outdir+fname.split('/')[-1]+".twas",'w')
    header=data[0].split()
    snp_index=header.index('SNP')
    if keep_columns==None:
        use_indices=range(len(header))
    if keep_columns!=None:
        use_indices=[header.index(i) for i in keep_columns if i in header]                
    a2_only=False
    if 'A1' in header:
        a2_only=True
        a1_index=header.index('A1')
        outf.write('\t'.join([header[i] for i in use_indices])+'\tA2\n') 
    else: 
        outf.write('\t'.join([header[i] for i in use_indices])+'\tA1\tA2\n') 
    for line in data[1::]:
        try:
            tokens=line.split()
            cur_snp=tokens[snp_index]
            if cur_snp in exclude_dict:
                continue 
            cur_freqs=freq_dict[cur_snp]  
            if (a2_only==True):
                if cur_freqs[0]==tokens[a1_index]:
                    a2=cur_freqs[1]
                else:
                    a2=cur_freqs[0]
                outf.write('\t'.join([tokens[i] for i in use_indices])+'\t'+a2+'\n')
            else:
                outf.write('\t'.join([tokens[i] for i in use_indices])+'\t'+'\t'.join(cur_freqs)+'\n')        
        except:
            print("WARNING! Skipping line:"+line)

def get_exclude_list(exclude_list_file):
    print("reading exclude SNP list") 
    exclude_data= open(exclude_list_file,'r').read().strip().split('\n')
    print("building exclusion dictionary") 
    exclude_dict=dict()
    for line in exclude_data:
        exclude_dict[line]=1
    del exclude_data 
    return exclude_dict

def main():
    args=parse_args()
    #read in the frequency information
    freq_dict=read_freqs(args.frequencies)
    file_list=open(args.gwas_output,'r').read().strip().split('\n')
    if args.exclude_list!=None: 
        exclude_dict=get_exclude_list(args.exclude_list)
    else:
        exclude_dict=dict() #nothing should be excluded unless specified as a command line argument
    outdir=args.outdir 
    if outdir.endswith('/')==False:
        outdir=outdir+'/' 
    for fname in file_list:
        annotate(fname,freq_dict,args.keep_columns,exclude_dict,outdir)
        
    
    

if __name__=="__main__":
    main()
    
