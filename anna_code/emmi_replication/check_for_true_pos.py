# identify the GWAS results from Ashley lab pipeline for the SNP's that Emmi's pipeline reports as significant. 
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="identify the GWAS results from Ashley lab pipeline for the SNP's that Emmi's pipeline reports as significant.")
    parser.add_argument("--gs_file")
    parser.add_argument("--gs_chrom",type=int)
    parser.add_argument("--gs_pos",type=int)
    parser.add_argument("--ashley_lab_gwas_prefix")
    parser.add_argument("--ashley_lab_chrom",type=int)
    parser.add_argument("--ashley_lab_pos",type=int)
    parser.add_argument("--outf")
    return parser.parse_args()
def main():
    args=parse_args()
    gs=open(args.gs_file,'r').read().strip().split('\n')
    gs_dict=dict()
    for line in gs[2::]:
        tokens=line.split()
        entry=tuple([tokens[args.gs_chrom],tokens[args.gs_pos]])
        gs_dict[entry]=line
    print("built dictionary of Emmi's hits")
    outf=open(args.outf,'w')
    header=["SNP","chr","pos","gene","effect_allele","other_allele","eaf","beta","se","p_value","beta","se","p_value","CHR","SNP","BP","A1","TEST","NMISS","BETA","STAT","P"]
    outf.write("\t".join(header)+"\n")
    for chrom in range(1,23):
        cur_file=open(args.ashley_lab_gwas_prefix+"."+str(chrom)+".continuous.assoc.linear",'r').read().strip().split('\n')
        for line in cur_file[1::]:
            tokens=line.split()
            try:
                chrom=tokens[args.ashley_lab_chrom] 
                pos=tokens[args.ashley_lab_pos]
                entry=tuple([chrom,pos])
                if entry in gs_dict:
                    outf.write(gs_dict[entry]+'\t'+'\t'.join(tokens)+'\n')
            except:
                print(str(tokens))
                


if __name__=="__main__":
    main()
    
