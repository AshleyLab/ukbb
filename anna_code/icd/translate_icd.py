
import argparse
import pdb
def parse_args():
    parser=argparse.ArgumentParser(description="tranlsate ICD code tallies")
    parser.add_argument("--icd_hist",default="ICD_tallies.tsv")
    parser.add_argument("--translation",default="ICDcoding.csv")
    parser.add_argument("--outf",default="ICD_tallies.translated.tsv")
    parser.add_argument("--index",default=0,type=int) 
    return parser.parse_args()

def main():
    args=parse_args() 
    translation=open(args.translation,'r').read().strip().split('\n')
    mapping=dict()
    for line in translation:
        tokens=line.split('\t')
        mapping[tokens[0]]=tokens[1]
    tallies=open(args.icd_hist,'r').read().strip().split('\n')
    outf=open(args.outf,'w')
    outf.write(tallies[0]+'\n') 
    for line in tallies[1::]:
        tokens=line.split('\t')
        try:
            mapped=mapping[tokens[args.index]]
            outf.write(line+'\t'+mapped+'\n')
        except:
            continue 
if __name__=="__main__":
    main()
    
