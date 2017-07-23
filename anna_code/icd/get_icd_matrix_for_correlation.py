#generates a subject x ICD code binary matrix for regression against phenotype features
import argparse
#import pdb 
def parse_args():
    parser=argparse.ArgumentParser(description="generates a subject x ICD code binary matrix for regression against phenotype features")
    parser.add_argument("--icd_list",default="codes_over_1000.txt")
    parser.add_argument("--ukbb_bulk_data",default="/oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb9419.tab")
    parser.add_argument("--outf",default="icd_matrix_for_regression.txt")
    return parser.parse_args()

def main():
    args=parse_args()
    codes=open(args.icd_list,'r').read().strip().split('\n')
    code_dict=dict()
    for code in codes:
        code_dict[code]=1

    dataf=open(args.ukbb_bulk_data,'r')
    line=dataf.readline()
    line=dataf.readline()
    counter=1
    outf=open(args.outf,'w')
    outf.write("IID"+'\t'+'\t'.join(codes)+'\n')

    while line!="":
        tokens=line.split('\t')
        subject=tokens[0]
        mat_for_reg=dict() 
        for code in code_dict:
            mat_for_reg[code]='0' 
        cur_codes=[i.replace('"','') for i in set(tokens[7:-2])]       
        for code in cur_codes:
            if code in code_dict:
                mat_for_reg[code]='1'
        outf.write(subject)
        for code in codes:
            outf.write('\t'+mat_for_reg[code])
        outf.write('\n') 
        line=dataf.readline()
        counter+=1
        if counter%1000==0:
            print(str(counter))
            
        
if __name__=="__main__":
    main()
    
