import argparse
def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument('--motif_hits')
    parser.add_argument('--outf')
    return parser.parse_args()

    
def main():
    args=parse_args()
    data=open(args.motif_hits,'r').read().strip().split('\n')
    outf=open(args.outf,'w')
    outf.write(data[0]+'\n')
    snp_dict=dict() 
    for line in data[1::]:
        tokens=line.split('\t')
        snp=tokens[2]
        motif_fam=tokens[3].split('/')[0].split('_')[1]
        motif=tokens[3].split('_')[1] 
        score=float(tokens[5])
        if snp not in snp_dict:
            snp_dict[snp]=dict()
        if motif_fam not in snp_dict:
            snp_dict[snp][motif_fam]=[score,line]
        else:
            existing_score=snp_dict[snp][motif_fam][0]
            if score > existing_score:
                snp_dict[snp][motif_fam]=[score,line]
    for snp in snp_dict:
        for motif_fam in snp_dict[snp]:
            outf.write(snp_dict[snp][motif_fam][1]+'\n')
            
    
if __name__=="__main__":
    main()
    
    
