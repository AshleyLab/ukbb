#converts data encoding from one UKBB format to another
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="")
    parser.add_argument("--data_file")
    parser.add_argument("--columns",type=int,nargs="+",help="0-indexed column with desired field of interest")
    parser.add_argument("--encoding_map",nargs="+")
    parser.add_argument("--outf")
    return parser.parse_args()

def make_conversion_dict(encoding_map):
    conversion_dict=dict()
    encoding=open(encoding_map,'r').read().strip().split('\n')
    for line in encoding:
        tokens=line.split('\t')
        conversion_dict[tokens[0]]=tokens[1] 
    return conversion_dict

def main():
    args=parse_args()
    conversion_dicts=[make_conversion_dict(i) for i in args.encoding_map]
    conversion_metadict=dict()
    for i in range(len(args.columns)):
        conversion_metadict[args.columns[i]]=conversion_dicts[i]        
    data=open(args.data_file,'r').read().strip().split('\n')
    outf=open(args.outf,'w')
    outf.write(data[0]+'\n')
    for line in data[1::]:
        tokens=line.split('\t')
        for index in args.columns:            
            toconvert=tokens[index]
            tokens[index]=conversion_metadict[index][toconvert]
        outf.write('\t'.join(tokens)+'\n')
        
if __name__=="__main__":
    main()
    
