#no disease=1, has disease=2
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="collapse icd matrix")
    parser.add_argument("--icd_matrix",default="icd_matrix_cardiovascular.txt")
    parser.add_argument("--outf",default="icd_phenotypes_binary.txt")
    return parser.parse_args()


def main():
    args=parse_args()
    data=open(args.icd_matrix,'r').read().strip().split('\n')
    print("loaded data matrix")
    outf=open(args.outf,'w')
    outf.write("FID\tIID\tCardiovascularDisease\n")
    counter=0
    for line in data[1::]:
        counter+=1
        if counter % 1000==0:
            print(str(counter))
        tokens=line.split('\t')
        iid=tokens[0]
        max_val=max([int(i) for i in tokens[1::]])
        if max_val==1:
            outf.write(iid+'\t'+iid+'\t'+'2'+'\n')
        else:
            outf.write(iid+'\t'+iid+'\t'+'1'+'\n')
            

if __name__=="__main__":
    main()
    
    
    
