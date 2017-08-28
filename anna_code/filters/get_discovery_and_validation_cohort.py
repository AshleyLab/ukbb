import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="split subjects into discovery & replication cohort")
    parser.add_argument("--subject_list")
    parser.add_argument("--percent_discovery",default=0.66,type=float)
    parser.add_argument("--out_prefix")
    return parser.parse_args()

def main():
    args=parse_args()
    subjects=open(args.subject_list,'r').read().strip().split('\n')
    discovery=[]
    validation=[]
    import random
    for i in range(len(subjects)):
        cur_line=subjects[i]
        new_random=random.random()
        if new_random < args.percent_discovery:
            discovery.append(cur_line)
        else:
            validation.append(cur_line)

    outf_discovery=open(args.out_prefix+".discovery",'w')
    outf_discovery.write('\n'.join(discovery))
    outf_validation=open(args.out_prefix+".validation",'w')
    outf_validation.write('\n'.join(validation))
    
    
if __name__=="__main__":
    main()
    
