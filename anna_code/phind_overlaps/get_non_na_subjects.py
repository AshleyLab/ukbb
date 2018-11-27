import argparse 
def parse_args(): 
    parser=argparse.ArgumentParser(description="identify subjects that have a non-NA value for any of the trait-associated fields in the file")
    parser.add_argument("--inputf") 
    parser.add_argument("--outf") 
    return parser.parse_args() 
def main() : 
    args=parse_args() 
    data=open(args.inputf,'r').read().strip().split('\n') 
    outf=open(args.outf,'w') 
    for line in data[1::]: 
        tokens=line.split('\t') 
        values=list(set(tokens[1::]))
        if len(values)>1: 
            outf.write(tokens[0]+'\n')
        elif values[0]!="NA": 
            outf.write(tokens[0]+'\n')

if __name__=="__main__": 
    main() 
