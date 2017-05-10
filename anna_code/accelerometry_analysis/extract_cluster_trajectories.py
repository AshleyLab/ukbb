import argparse
from Params import *
import pdb 
def parse_args():
    parser=argparse.ArgumentParser("get mean trajectory for subjects in specified cluster")
    parser.add_argument("--base_dir")
    parser.add_argument("--subject_to_cluster")
    parser.add_argument("--cluster_of_interest",default=None)
    parser.add_argument("--outf")
    return parser.parse_args()

def load_map(mapfile):
    data=open(mapfile,'r').read().strip().split('\n')
    subject_to_cluster=dict()
    header=data[0].split('\t')
    print(str(header))
    subject_index=header.index('subject')
    group_index=header.index('group')    
    for line in data[1::]:
        tokens=line.split('\t')
        cur_subject=tokens[subject_index]
        cur_cluster=tokens[group_index]
        if cur_cluster not in subject_to_cluster:
            subject_to_cluster[cur_cluster]=[cur_subject]
        else:
            subject_to_cluster[cur_cluster].append(cur_subject)
    return subject_to_cluster

def main():
    args=parse_args()
    subject_to_cluster=load_map(args.subject_to_cluster)
    clusters_of_interest=subject_to_cluster.keys()
    if args.cluster_of_interest!=None:
        clusters_of_interest=[args.cluster_of_interest]
    outf=open(args.outf,'w')
    for cluster in clusters_of_interest:
        print("pooling for cluster:"+str(cluster))
        subjects=subject_to_cluster[cluster]
        num_entries=len(subjects)
        print(str(num_entries))
        total_vals=[0]*expected_number_datapoints
        num_processed=0
        for subject in subjects:
            num_processed+=1
            #print(str(subject))
            if num_processed %10==0:
                print(str(num_processed))
            subject_subdir='/'+subject[0]+'/'
            subject_data=open(args.base_dir+subject_subdir+subject+"_90004_0_0.csv",'r').read().strip().split('\n')
            #handle missing values -- replace with 0 if no measurement ? TODO: is this the best way to handlel missing values? 
            subject_data=[0 if x=="" else float(x) for x in subject_data]
            #increment the toal 
            total_vals=[x + y for x, y in zip(total_vals,subject_data)]
        #divide by the number of entries to get the mean trajectory
        total_vals=[i/num_entries for i in total_vals]
        outf.write(cluster+'\t'+str(num_entries)+'\t'+'\t'.join([str(i) for i in total_vals])+'\n')

    
if __name__=="__main__":
    main()
    
