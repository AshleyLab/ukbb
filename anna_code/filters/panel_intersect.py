#check to see the overlap of the directly genotyped (called) panel with the imputed panel.
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="check overlap between imputed and directly genotyped panels")
    parser.add_argument("--imputed_snps")
    parser.add_argument("--called_snps")
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    called=open(args.called_snps,'r').read().strip().split('\n')
    called_dict=dict()
    for line in called:
        tokens=line.split('\t')
        pos=tokens[3]
        called_dict[pos]=line 
    print("read in called snps:"+str(len(called_dict.keys())))
    outf=open(args.outf,'w')
    imputed=open(args.imputed_snps,'r')
    cur_imputed=imputed.readline()
    while cur_imputed!="":
        tokens=cur_imputed.split('\t')
        pos=tokens[1]
        if pos in called_dict:
            called_dict.__delitem__(pos)
        cur_imputed=imputed.readline()
    print("called snps not on imputed panel:"+str(len(called_dict.keys())))
    outf.write('\n'.join(called_dict.values()))
    
    
    
if __name__=="__main__":
    main()
    
