#tallies up the ICD codes and maps them to the corresponding disease meaning
import argparse
import pdb
def parse_args():
    parser=argparse.ArgumentParser(description="tallies up the ICD10 codes in UKBB subjects and maps them to corresponding disease meaning")
    parser.add_argument("--bulk_data",default="/oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb9419.tab")
    parser.add_argument("--translation",default="ICDcoding.csv")
    parser.add_argument("--outf",default="ICD_tallies.tsv")
    return parser.parse_args()

def main():
    args=parse_args()
    code_hist=dict()
    dataf=open(args.bulk_data,'r')
    line=dataf.readline()#we skip the first line because it's the header
    line=dataf.readline()
    counter=1
    while line!="":
        tokens=line.split('\t')
        codes=set(tokens[7:-2])
        for code in codes:
            if code not in code_hist:
                code_hist[code]=1
            else:
                code_hist[code]+=1
        line=dataf.readline()
        counter+=1
        if counter%1000==0:
            print(str(counter))
    outf=open(args.outf,'w')
    for code in code_hist:
        outf.write(code+'\t'+str(code_hist[code])+'\n')
        
    
    
        
if __name__=="__main__":
    main()
    
    
