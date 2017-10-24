#filter and combine the clumped results.
#filter by panel
#filter by p-value
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="filter outputs of PLINK clump by panel and p-value")
    parser.add_argument("--sourcef")
    parser.add_argument("--pval_thresh",type=float,default=5e-8)
    parser.add_argument("--bad_snp_list",default="/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/hrc/bad.snps.txt")
    parser.add_argument("--outf")
    parser.add_argument("--verify_parts",nargs="+")
    return parser.parse_args()

def make_bad_snp_dict(bad_snp_list):
    data=open(bad_snp_list,'r').read().strip().split('\n')
    bad_snp_dict=dict()
    for line in data:
        bad_snp_dict[line]=1
    print("made bad snp dict")
    return bad_snp_dict 

def get_sig_snps(parts,thresh):
    sig_snp_dict=dict()
    for chrom in range(1,23):
        data=open(parts[0]+str(chrom)+parts[1],'r').read().strip().split('\n')
        for line in data[1::]:
            tokens=line.split( )
            pval=tokens[-1]
            if pval!="NA":
                pval=float(pval)
                if pval<thresh:
                    snp=tokens[2]
                    sig_snp_dict[snp]=pval
    print("generated dictionary of significant SNPs") 
    return sig_snp_dict 

def filter_clump(clump,bad_snps,sig_snps):
    min_p=1
    index_snp=None
    for snp in clump:
        if snp in bad_snps:
            continue
        if snp in sig_snps:
            cur_pval=sig_snps[snp]
            if cur_pval<min_p:
                min_p=cur_pval
                index_snp=snp
    return index_snp,min_p 

def main():
    args=parse_args()
    bad_snps=make_bad_snp_dict(args.bad_snp_list)
    sig_snps=get_sig_snps(args.verify_parts,args.pval_thresh) 
    data=open(args.sourcef,'r').read().strip().split('\n')
    outf=open(args.outf,'w')
    outf.write(data[0])
    for line in data:
        if line.startswith("CHR"):
            continue 
        tokens=line.split()
        try:
            clump=[tokens[2]]+[i.split('(')[0] for i in tokens[-1].split(',')]
        except:
            print(str(tokens))
            continue
        print(str(clump)) 
        index_snp,index_snp_p=filter_clump(clump,bad_snps,sig_snps)
        if index_snp!=None:
            outf.write(index_snp+'\t'+str(index_snp_p)+'\n')
if __name__=="__main__":
    main()
    
