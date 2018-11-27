import argparse 
def parse_args(): 
    parser=argparse.ArgumentParser("Get intersection of subject files") 
    parser.add_argument("--files_to_intersect",nargs="+") 
    parser.add_argument("--outf") 
    return parser.parse_args() 

def main(): 
    args=parse_args() 
    cur_intersection=set(open(args.files_to_intersect[0]).read().strip().split('\n'))
    for i in range(1,len(args.files_to_intersect)): 
        data=set(open(args.files_to_intersect[i]).read().strip().split('\n'))
        cur_intersection=data.intersection(cur_intersection)
    outf=open(args.outf,'w') 
    cur_intersection='\n'.join(list(cur_intersection))
    outf.write(cur_intersection+'\n')


if __name__=="__main__": 
    main() 

