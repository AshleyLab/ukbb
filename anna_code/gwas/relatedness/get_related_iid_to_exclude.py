#finds all iid's with kinship over a threshold and randomly picks one individual to exclude.
import argparse
import random 
def parse_args():
    parser=argparse.ArgumentParser(description="finds all iid's with kinship over a threshold and randomly picks one individual to exclude.")
    parser.add_argument("--relationship_info",default="/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/ukb2228_rel_s488374.dat")
    parser.add_argument("--kinship_thresh",default=0.22,type=float)
    parser.add_argument("--outf",default="kinship_higher_0.22.txt")
    return parser.parse_args()

def main():
    args=parse_args()
    rel=open(args.relationship_info,'r').read().strip().split('\n')
    outf=open(args.outf,'w')
    removed=dict() 
    thresh=args.kinship_thresh
    for line in rel[1::]:
        tokens=line.split()
        kinship=float(tokens[-1])
        if kinship >=thresh:
            #check if one subject in the pair has already been removed
            subject1=tokens[0]
            subject2=tokens[1]
            if subject1 in removed:
                continue
            elif subject2 in removed:
                continue
            else: 
                coin_flip=random.random()
                if coin_flip <0.5: 
                    outf.write(subject1+'\n')
                    removed[subject1]=1
                else:
                    outf.write(subject2+'\n')
                    removed[subject2]=1

if __name__=="__main__":
    main()
    
