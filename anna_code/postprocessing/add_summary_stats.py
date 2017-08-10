#add number of traits, mean effect size, mean score
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
    outf.write('SNP\tNumberOfFeaturesWithSignificantAssoc\tMeanZscore\tMeanEffectSize\tFeatures\n')
    for line in data[1::]:
        tokens=line.split('\t')
        snp=tokens[1]
        features_affected=[] 
        effect=[tokens[i] for i in range(10,len(tokens),3)]
        num_p_vals=0 
        for i in range(10,len(tokens),3):
            if tokens[i]!="":
                features_affected.append(header[i])
                num_p_vals+=1
                
        stat=[tokens[i] for i in range(11,len(tokens),3)] 
        while "" in effect:
            effect.remove("")
        while "" in stat:
            stat.remove("")
        effect=[float(i) for i in effect]
        stat=[float(i) for i in stat]
        effect=sum(effect)/len(effect)
        stat=sum(stat)/len(stat)
        outf.write(snp+'\t'+str(num_p_vals)+'\t'+str(stat)+'\t'+str(effect)+'\t'+','.join(features_affected)+'\n')
        
if __name__=="__main__":
    main()
    
