#generate a matrix for analysis with R upsetR package 
import argparse 
def parse_args(): 
    parser=argparse.ArgumentParser(description="") 
    parser.add_argument("--lists",nargs="+") 
    parser.add_argument("--labels",nargs="+")
    parser.add_argument("--outf") 
    return parser.parse_args() 

def main(): 
    args=parse_args() 
    outf=open(args.outf,'w') 
    subject_dict=dict() 
    all_subjects=set([]) 
    for i in range(len(args.lists)): 
        cur_label=args.labels[i] 
        cur_list=open(args.lists[i]).read().strip().split('\n')  
        all_subjects=all_subjects.union(set(cur_list))
        subject_dict[cur_label]=dict() 
        for entry in cur_list: 
            subject_dict[cur_label][entry]=1 
        print(cur_label+":done")
    labels=list(subject_dict.keys())
    outf.write("Subject"+'\t'+'\t'.join(labels)+'\n')
    all_subjects=list(all_subjects) 
    for subject in all_subjects: 
        outf.write(subject) 
        for label in labels: 
            if subject in subject_dict[label]: 
                outf.write('\t1') 
            else: 
                outf.write('\t0') 
        outf.write('\n')

if __name__=="__main__": 
    main() 
