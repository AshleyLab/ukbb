import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="check if extra SNPs in our pipeline are in LD w/ called SNPs")
    parser.add_argument("--hp_dir")
    parser.add_argument("--gs_snps")
    parser.add_argument("--called_snps")
    return parser.parse_args()

def main():
    args=parse_args()
    gs_snps=open(args.gs_snps,'r').read().strip().split('\n')
    gs_dict=dict()
    for snp in gs_snps:
        gs_dict[snp]=1 
    called_snps=set(open(args.called_snps,'r').read().strip().split('\n'))
    gs_snps=set(gs_snps)
    
    in_same_block=set([]) 
    for chrom in range(1,23):
        hp=open(args.hp_dir+'/'+'plink.blocks.'+str(chrom)+'.blocks','r').read().strip().split('\n')
        for line in hp:
            tokens=line.split()[1::]
            for t in tokens:
                if t in gs_dict:
                    #add to set of snps in same block
                    tokens=set(tokens)
                    in_same_block=in_same_block.union(tokens)
                    break
    print('made same block dict')
    unknown_snps=called_snps-gs_snps
    matched=unknown_snps.intersection(in_same_block)
    unmatched=unknown_snps - in_same_block
    print("IN LD WITH EMMI'S SNPS:"+str(len(matched)))
    print("NOT IN LD WITH EMMI'S SNPS:"+str(len(unmatched)))
    
if __name__=="__main__":
    main()
    
