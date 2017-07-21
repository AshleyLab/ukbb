#Removes David H's excluded subjects from the pca-filtered subject list
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="removes subjects on exclusion list from pca-filtered subject list")
    parser.add_argument("--pca_filtered_list")
    parser.add_argument("--exclusion_list")
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    exclusion_list=open(args.exclusion_list,'r').read().strip().split('\n')
    exclusion_dict=dict()
    for line in exclusion_list:
        tokens=line.split()
        exclusion_dict[tokens[0]]=1
    pca_filtered_list=open(args.pca_filtered_list,'r').read().strip().split('\n')
    outf=open(args.outf,'w')
    for line in pca_filtered_list:
        tokens=line.split()
        if tokens[0] not in exclusion_dict:
            outf.write(line+'\n')
            
if __name__=="__main__":
    main()
    
