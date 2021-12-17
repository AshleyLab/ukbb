#add fields to a covariates or outcomes file by merging on FID from another file. 
import argparse 
def parse_args(): 
    parser=argparse.ArgumentParser(description="add fields to a covariates or outcomes file by merging on FID from another file. ")
    parser.add_argument("--mainf_f") 
    parser.add_argument("--source_f",nargs="+") 
    parser.add_argument("--source_fields",nargs="+") 
    parser.add_argument("--outf") 
    return parser.parse_args() 

def main(): 
    args=parse_args() 
    #merge by subjects 
    outf=open(args.outf,'w') 
    num_fields=len(args.source_fields)
    new_field_dict=dict() 
    for i in range(num_fields): 
        data=open(args.source_f[i],'r').read().strip().split('\n') 
        curfield=args.source_fields[i] 
        header=data[0].split('\t') 
        field_index=header.index(curfield) 
        new_field_dict[field_index]=dict() 
        for line in data[1::]: 
            tokens=line.split('\t') 
            subject=tokens[0] 
            fieldval=tokens[field_index] 
            new_field_dict[field_index][subject]=fieldval
    data=open(args.mainf_f,'r').read().strip().split('\n') 
    header=data[0] 
    outf.write(header+'\t'.join(args.source_fields)+'\n')
    for line in data[1::]: 
        tokens=line.split('\t') 
        subject=tokens[0] 
        outf.write(line) 


if __name__=="__main__": 
    main() 

