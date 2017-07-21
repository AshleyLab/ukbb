import pandas as pd
import argparse
import pdb

def parse_args():
    parser=argparse.ArgumentParser(description="collapses pairs of highly correlated features by selecting one feature from the pair to remove")
    parser.add_argument("--corr_mat")
    parser.add_argument("--thresh",type=float)
    parser.add_argument("--out_prefix")
    parser.add_argument("--must_keep",nargs="*")
    return parser.parse_args()

def main():
    args=parse_args()
    data=pd.DataFrame.from_csv(args.corr_mat,sep='\t',header=0)
    thresh=args.thresh
    remove=dict()
    keep=dict()
    must_keep=args.must_keep
    must_keep_dict=dict()
    for entry in must_keep:
        must_keep_dict[entry]=1
    print(str(must_keep_dict))
    features=data.index
    for f1 in features:
        for f2 in features:
            if f1==f2:
                continue
            if abs(data[f1][f2])>thresh:
                #remove one of the features!
                if ((f1 in remove) or( f2 in remove)):
                    #one of the features has already been removed, no need to do anything else 
                    continue
                elif (f1 not in must_keep_dict): 
                    #remove one of the features if it's not in our 'must keep' list 
                    remove[f1]=1
                elif (f2 not in must_keep_dict):
                    remove[f2]=1
                    
    for feature in features:
        if feature not in remove:
            keep[feature]=1
    data_subset=data.loc[keep.keys()][keep.keys()]
    data_subset.to_csv(args.out_prefix,sep="\t")
    outf=open(args.out_prefix+".header",'w')
    outf.write('\n'.join(keep.keys()))
    
    
if __name__=="__main__":
    main()
    
