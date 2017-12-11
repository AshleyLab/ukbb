import sys
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="consolidate phenotypes for MR analysis")
    parser.add_argument("--source_pheno",nargs="+")
    parser.add_argument("--features")
    parser.add_argument("--outf")
    parser.add_argument("--subjects") 
    return parser.parse_args()

def main():
    args=parse_args()
    features=open(args.features,'r').read().strip().split('\n')
    subjects=open(args.subjects,'r').read().strip().split('\n')
    
    feature_dict=dict()
    for feature in features:
        feature_dict[feature]=1
        
    subject_to_pheno=dict()
    for ph in args.source_pheno:
        data=open(ph,'r').read().strip().split('\n')
        header=data[0].split('\t')
        tokeep=[]
        for i in range(len(header)):
            cur_feature=header[i]
            if cur_feature in feature_dict:
                tokeep.append(i)
        for line in data[1::]:
            tokens=line.split('\t')
            subject=tokens[0]
            for index in tokeep:
                cur_feat=header[index]
                cur_val=tokens[index]
                if cur_val=="-1000":
                    cur_val="NA" 
                if cur_feat not in subject_to_pheno:
                    subject_to_pheno[cur_feat]=dict()
                subject_to_pheno[cur_feat][subject]=cur_val
    outf=open(args.outf,'w')
    outf.write('Subject')
    for feature in features:
        outf.write('\t'+feature)
    outf.write('\n') 
    for subject in subjects:
        outf.write(subject)
        for feature in features:
            if subject in subject_to_pheno[feature]:
                outf.write('\t'+subject_to_pheno[feature][subject])
            else:
                outf.write('\tNA') 
        outf.write('\n')
        
if __name__=="__main__":
    main()
    
