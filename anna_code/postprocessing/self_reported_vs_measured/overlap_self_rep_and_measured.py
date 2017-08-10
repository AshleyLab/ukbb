import argparse
import pdb
def parse_args():
    parser=argparse.ArgumentParser(description="which snp's are associated with both self-reported and measured values?")
    parser.add_argument("--self_reported_traits")
    parser.add_argument("--measured_traits")
    parser.add_argument("--top_snps")
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    self_reported=open(args.self_reported_traits,'r').read().strip().split('\n')
    self_rep_dict=dict()
    for entry in self_reported:
        self_rep_dict[entry]=1
    measured=open(args.measured_traits,'r').read().strip().split('\n')
    measured_dict=dict()
    for entry in measured:
        measured_dict[entry]=1
    data=open(args.top_snps,'r').read().strip().split('\n')
    header=data[0].split('\t')
    feat_index=header.index('Features') 
    outf=open(args.outf,'w')
    outf.write(data[0]+'\tMeasurementType'+'\n')
    for line in data[1::]:
        tokens=line.split('\t')
        print(str(tokens))
        assoc_features=tokens[feat_index].split(',')
        self_rep=False
        measured=False
        for f in assoc_features:
            if f in self_rep_dict:
                self_rep=True
            if f in measured_dict:
                measured=True
        if ((self_rep==True) and (measured==True)):
            outf.write(line+'\t'+'BOTH'+'\n')
        elif(self_rep==True):
            outf.write(line+'\t'+'Self-reported'+'\n')
        elif(measured==True):
            outf.write(line+'\t'+'Measured'+'\n')
        else:
            print(str(assoc_features))
            pdb.set_trace() 
            raise Exception()
            

if __name__=="__main__":
    main()
    
