#aggregate components of the mean trajectory
import argparse
from Params import *

def parse_args():
    parser=argparse.ArgumentParser("aggregate components of the mean trajectory")
    parser.add_argument("--base_name")
    parser.add_argument("--start_index",type=int)
    parser.add_argument("--end_index",type=int)
    parser.add_argument("--increment",type=int) 
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    start_index=args.start_index
    end_index=args.end_index
    increment=args.increment 
    total_vals=0
    summed_vals=[0]*expected_number_datapoints 
    while start_index < end_index:
        print(str(start_index)+":"+str(start_index+increment))
        try:
            data=open(args.base_name+'.'+str(start_index)+'.'+str(start_index+increment),'r').read().strip().split('\n')[0]
        except:
            increment=1000
            data=open(args.base_name+'.'+str(start_index)+'.'+str(start_index+increment),'r').read().strip().split('\n')[0]
        data=data.split('\t') 
        cluster_name=int(data[0])
        total_vals+=int(data[1])
        values=[float(i) for i in data[2::]]
        assert(len(values)==len(summed_vals)) #make sure the lengths match!
        summed_vals=[x + y for x, y in zip(summed_vals,values)]
        start_index+=increment
        
    #divide by the total number of entries to get the mean value
    summed_vals=[float(i)/total_vals for i in summed_vals]
    outf=open(args.outf,'w')
    outf.write(str(cluster_name)+'\t'+str(total_vals)+'\t'+'\t'.join([str(i) for i in summed_vals])+'\n')
    
if __name__=="__main__":
    main()
    
