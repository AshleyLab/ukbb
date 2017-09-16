bad_snps=open("/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/hrc/bad.snps.txt").read().strip().split('\n') 
bad_dict=dict()
for line in bad_snps:
    bad_dict[line]=1
del bad_snps 
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="filter SNPs on the panel to exclude")
    parser.add_argument("--input")
    return parser.parse_args()

def main():
    args=parse_args()
    outf=open(args.input+".filtered",'w')
    data=open(args.input,'r').read().strip().split('\n')
    for line in data:
        tokens=line.split()
        rs=tokens[1]
        if rs not in bad_dict:
            outf.write(line+'\n')
            

if __name__=="__main__":
    main()
    
    
    
