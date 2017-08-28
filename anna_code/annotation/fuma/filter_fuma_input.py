#filter fuma input file by maf or info columns
import argparse
import pdb
def parse_args():
    parser=argparse.ArgumentParser(description="filter fuma inputs by maf and/or info")
    parser.add_argument("--fuma_input")
    parser.add_argument("--maf_filter_upper",default=None,type=float)
    parser.add_argument("--maf_filter_lower",default=None,type=float)
    parser.add_argument("--info_filter",default=None,type=float)
    parser.add_argument("--metadata")
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    metadata=open(args.metadata,'r').read().strip().split('\n')
    header=metadata[0].split('\t')
    maf_index=header.index('MAF')
    info_index=header.index('Info')
    snp_index=header.index('SNP') 
    metadata_dict=dict()
    for line in metadata[1::]:
        tokens=line.split('\t')
        #print(str(tokens))
        snp=tokens[snp_index]
        maf=float(tokens[maf_index])
        info=float(tokens[info_index])
        metadata_dict[snp]=[maf,info]
    print("made top snp reference dictionary")
    outf=open(args.outf,'w')
    fuma_input=open(args.fuma_input,'r').read().strip().split('\n')
    outf.write(fuma_input[0]+'\n')
    for line in fuma_input[1::]:
        tokens=line.split('\t')
        snp=tokens[1]
        anno=metadata_dict[snp]
        cur_maf=anno[0]
        cur_info=anno[1]
        if ((args.maf_filter_upper!=None) and (args.info_filter!=None)):
            if (cur_maf >=args.maf_filter_lower):
                if (cur_maf <=args.maf_filter_upper):
                    if (cur_info >=args.info_filter):
                        outf.write(line+'\n')
        elif (args.maf_filter_upper!=None):
            #pdb.set_trace() 
            if (cur_maf >=args.maf_filter_lower):
                if (cur_maf <=args.maf_filter_upper): 
                    outf.write(line+'\n')
        elif (args.info_filter!=None):
            if (cur_info >= args.info_filter):
                outf.write(line+'\n')
        else:
            outf.write(line+'\n')
        
    

if __name__=="__main__":
    main()
    
    
