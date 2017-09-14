import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="generate flank bed for SNP hits")
    parser.add_argument("--snp_hit_file")
    parser.add_argument("--outf")
    parser.add_argument("--flank_size",type=int) 
    return parser.parse_args()

def make_flank_bed(data,flank_size,outf):
    outf_discovery=open(outf+".discovery.bed",'w')
    outf_replication=open(outf+".replication.bed",'w')
    for line in data:
        tokens=line.split('\t')
        #print(str(tokens))
        discovery_chrom=tokens[0]
        if discovery_chrom!="":
            d_rs=tokens[1]
            d_pos=int(tokens[2])
            d_startpos=d_pos-flank_size
            if d_startpos <1:
                d_startpos=1
            d_endpos=d_pos+flank_size
            outf_discovery.write(discovery_chrom+'\t'+str(d_startpos)+'\t'+str(d_endpos)+'\t'+d_rs+'\n')
        if (len(tokens)<10):
            continue 
        replication_chrom=tokens[9]
        if replication_chrom!="":
            r_rs=tokens[10]
            r_pos=int(tokens[11])
            r_startpos=r_pos-flank_size
            if r_startpos <1:
                r_startpos=1
            r_endpos=r_pos+flank_size
            outf_replication.write(replication_chrom+'\t'+str(r_startpos)+'\t'+str(r_endpos)+'\t'+r_rs+'\n')

def main():
    args=parse_args()
    data=open(args.snp_hit_file,'r').read().strip().split('\n')
    make_flank_bed(data,args.flank_size,args.outf)
    
if __name__=="__main__":
    main()
    
