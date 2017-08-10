#gets number of significant variants (thresholded by maf) for each feature
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="add annotation stats to top snp summary")
    parser.add_argument("--top_snps")
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    data=open(args.top_snps,'r').read().strip().split('\n')
    outf=open(args.outf,'w')
    header=data[0].split('\t') 
    outf.write('Feature\tNumSNPs\tNumSNPsMafBelow0.01\n')
    feature_to_snp=dict()
    
    for line in data[1::]:
        tokens=line.split('\t')
        snp=tokens[1]
        try:
            maf=float(tokens[7])
        except:
            print(str(tokens))
        features_affected=[] 
        effect=[tokens[i] for i in range(14,len(tokens),3)]
        num_p_vals=0 
        for i in range(14,len(tokens),3):
            if tokens[i]!="":
                features_affected.append(header[i].split('.')[0])
                num_p_vals+=1
        for f in features_affected:
            if f not in feature_to_snp:
                feature_to_snp[f]=[0,0]
            if maf>=0.01:
                feature_to_snp[f][1]+=1
            feature_to_snp[f][0]+=1
    for feature in feature_to_snp:
        outf.write(feature+'\t'+str(feature_to_snp[feature][0])+'\t'+str(feature_to_snp[feature][1])+'\n')
        
if __name__=="__main__":
    main()
    
