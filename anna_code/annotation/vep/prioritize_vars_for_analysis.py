#ranks the VEP results by impact and maf for checking against the literature
import pdb 
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="Rank VEP results by impact and MAF")
    parser.add_argument("--vep_results")
    parser.add_argument("--top_snps")
    parser.add_argument("--impact_column",type=int)
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    top_snps=open(args.top_snps,'r').read().strip().split('\n')
    header=top_snps[0].split('\t')
    maf_index=header.index('MAF')
    chrom_index=header.index('CHR')
    pos_index=header.index('POS')
    top_snps_dict=dict()
    for line in top_snps[1::]:
        tokens=line.split('\t')
        chrom=tokens[chrom_index]
        pos=tokens[pos_index]
        maf=tokens[maf_index]
        cur_entry=tuple([chrom,pos])
        top_snps_dict[cur_entry]=float(maf)

    print("made frequency dictionary for top SNPs")
    vep_results=open(args.vep_results,'r').read().strip().split('\n')
    outf=open(args.outf,'w')
    header=vep_results[0].split('\t')
    outf.write('MAF'+'\t'+'\t'.join(header)+'\n')
    for line in vep_results[1::]:
        tokens=line.split('\t')
        snp_entry=tuple(tokens[0].split('_')[0:2])
        impact=tokens[args.impact_column]
        if impact in ["HIGH","MODERATE"]:
            curmaf=top_snps_dict[snp_entry]
            if curmaf >=0.01:
                #KEEP!
                outf.write(str(curmaf)+'\t'+line+'\n')
if __name__=="__main__":
    main()
    
    
