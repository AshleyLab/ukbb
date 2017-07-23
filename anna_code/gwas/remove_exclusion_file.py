#Removes David H's excluded subjects from the pca-filtered subject list
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="removes subjects on exclusion list from pca-filtered subject list")
    parser.add_argument("--pca_filtered_list")
    parser.add_argument("--exclusion_list",nargs="+")
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    pca_filtered_list=open(args.pca_filtered_list,'r').read().strip().split('\n')
    outf=open(args.outf,'w')
    exclusion_dict=dict()
    for l in args.exclusion_list: 
        exclusion_list=open(l,'r').read().strip().split('\n')
        for line in exclusion_list:
            tokens=line.split()
            exclusion_dict[tokens[0]]=1
    for line in pca_filtered_list:
        tokens=line.split()
        if tokens[0] not in exclusion_dict:
            outf.write(line+'\n')
            
if __name__=="__main__":
    main()
    
