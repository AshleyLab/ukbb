#generate file to download ukbb data in batches of 500
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="generate file to download ukbb data in batches of 500")
    parser.add_argument("--sourcefile")
    parser.add_argument("--numthreads",type=int)
    parser.add_argument("--increment",type=int,default=500) 
    return parser.parse_args()

def main():
    args=parse_args()
    data=open(args.sourcefile,'r').read().strip().split('\n')
    num_subjects=len(data)
    num_threads=args.numthreads
    increment=args.increment
    outf_files=[]
    for t in range(num_threads):
        cur_out=open('download_data.'+str(t)+'.sh','w')
        cur_out.write('#/bin/bash\n')
        outf_files.append(cur_out)
    print("set up output files!")
    cur_index=1
    out_index=0
    log_index=1
    while cur_index <num_subjects:
        outf_files[out_index].write('/srv/scratch/annashch/ukbb/ukbfetch -bto_download.accel.txt -s'+str(cur_index)+' -m'+str(increment)+' -of'+str(log_index)+'\n')
        cur_index+=increment
        out_index+=1
        if out_index>=len(outf_files):
            out_index=0
        log_index+=1

if __name__=="__main__":
    main()
