#compute intersection of availabel genotype & phenotype data, filter phenotype data frame accordingly
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="filter phenotype data to available genotype data")
    parser.add_argument('--input_df')
    parser.add_argument('--subject_ids')
    parser.add_argument('--outf')
    return parser.parse_args()

def main():
    args=parse_args()
    subject_dict=dict()
    for subject in open(args.subject_ids,'r').read().strip().split('\n'):
        subject_dict[subject]=1
    print("built subject dictionary")
    outf=open(args.outf,'w')
    data=open(args.input_df,'r').read().strip().split('\n')
    header=data[0]
    outf.write(header+'\n')
    for line in data[1::]: 
        tokens=line.split('\t')
        subject=tokens[0]
        if subject in subject_dict:
            outf.write(line+'\n')
if __name__=="__main__":
    main()
    
