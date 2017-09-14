#adds the effect size and associated feature information for the prioritized SNP hits.
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="add effect size and associated feature information for the prioritized SNP hits.")
    parser.add_argument("--snp_hit_file")
    parser.add_argument("--effect_col_name")
    parser.add_argument("--feature_col_name")
    parser.add_argument("--file_to_annotate") 
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    data=open(args.snp_hit_file,'r').read().strip().split('\n')
    header=data[0].split('\t')
    effect_col_index=header.index(args.effect_col_name)
    feature_col_index=header.index(args.feature_col_name)
    effect_dict=dict()
    for line in data[1::]:
        tokens=line.split('\t')
        chrom=tokens[0]
        pos=tokens[2]
        snp_entry=tuple([chrom,pos])
        cur_feature=tokens[feature_col_index]
        cur_effect=tokens[effect_col_index]
        if snp_entry not in effect_dict:
            effect_dict[snp_entry]=dict() 
        effect_dict[snp_entry][cur_feature]=cur_effect
    print("made effect dict")
    annotation=open(args.file_to_annotate,'r').read().strip().split('\n')
    header=annotation[0]
    outf=open(args.outf,'w')
    outf.write('Feature\tEffectSize\t'+header+'\n')
    for line in annotation[1::]:
        tokens=line.split('\t')
        snp_entry=tuple(tokens[1].split('_')[0:2])
        feature_hits=effect_dict[snp_entry]
        for feature in feature_hits:
            outf.write(feature+'\t'+feature_hits[feature]+'\t'+line+'\n')
            
    
if __name__=="__main__":
    main() 
