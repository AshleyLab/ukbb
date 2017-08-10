import argparse
import glob

def parse_args():
    parser=argparse.ArgumentParser(description="Find GWAS tasks that should be re-run")
    parser.add_argument("--base_dir")
    parser.add_argument("--task_list")
    parser.add_argument("--rerun_tasks")
    parser.add_argument("--rerun_chroms")
    return parser.parse_args()

def main():
    args=parse_args()
    outf_tasks=open(args.rerun_tasks,'w')
    outf_chroms=open(args.rerun_chroms,'w')
    tasks=open(args.task_list,'r').read().strip().split('\n')
    for task in tasks:
        for chrom in range(1,23):
            fname=args.base_dir+'/'+task+'/'+task+'.'+str(chrom)+'.'+'*adjusted'
            hits=glob.glob(fname)
            if len(hits)==0:
                outf_tasks.write(task+'\n')
                outf_chroms.write(str(chrom)+'\n')
                
    
if __name__=="__main__":
    main()
    
