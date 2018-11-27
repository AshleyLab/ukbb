import argparse 
def parse_args(): 
    parser=argparse.ArgumentParser(" Extract field values for subjects for specified list of fields")
    parser.add_argument("--table") 
    parser.add_argument("--fields",nargs="+") 
    parser.add_argument("--outf") 
    return parser.parse_args() 
def main(): 
    args=parse_args() 
    outf=open(args.outf,'w')
    indices=[] 
    #open file for reading 
    with open(args.table,'r') as fp: 
        #iterate through file line-by-line to keep memory footprint low 
        for cnt,line in enumerate(fp): 
            tokens=line.split() 
            if cnt % 1000==0: 
                print(cnt)
            if cnt==0: 
                #get the header indices
                for f in args.fields:
                    indices.append(tokens.index(f))
                print(indices)
                #write the output header 
                outf.write("Subject"+'\t'+'\t'.join(args.fields)+'\n')
            else: 
                #extract subject id and specified fields from the line
                subject=tokens[0] 
                #write the field values for the current subject 
                outf.write(subject+'\t'+'\t'.join(tokens[i] for i in indices)+'\n')

if __name__=="__main__": 
    main() 
