#get the number of transition states for each subject at various mg thresholds
import argparse
import pdb 
def parse_args():
    parser=argparse.ArgumentParser(description="get the number of transition states for each subject at various mg thresholds")
    parser.add_argument("--subject_file")
    parser.add_argument("--mg_cutoff",type=int)
    parser.add_argument("--outf")
    parser.add_argument("--consecutive_val_thresh",type=int,default=2)
    parser.add_argument("--start_index",type=int,default=None)
    parser.add_argument("--end_index",type=int,default=None) 
    return parser.parse_args()

def compress(data,mg_cutoff,cons_thresh):
    compressed_data=[]
    num_high=0
    num_low=0 
    for entry in data:
        try:
            value=float(entry.split(',')[0])
        except:
            #print(entry)
            continue
        if value >= mg_cutoff:
            if num_low==0: #we are in a run of high values 
                num_high+=1
            else: # We are not in a run of high values
                if num_low > cons_thresh: 
                    #add a record of low values
                    compressed_data.append(-1)
                #reset!
                num_high=1
                num_low=0
        else:
            if num_high==0: #We are in a run of low values
                num_low+=1
            else: #We are not in a run of low values
                if num_high > cons_thresh:
                    #add a record of high values
                    compressed_data.append(1)
                #reset!
                num_low=1
                num_high=0
    #record the final run of values!
    if num_high > cons_thresh:
        compressed_data.append(1)
    if num_low > cons_thresh:
        compressed_data.append(-1) 
    return compressed_data

def get_number_transition_states(data):
    num_transitions=0
    if len(data)==0:
        return num_transitions 
    last_elem=data[0]
    for elem in data[1::]:
        if elem!=last_elem:
            num_transitions+=1
        last_elem=elem 
    return num_transitions

def main():
    args=parse_args()
    outf=open(args.outf,'w')
    outf.write('Subject\tTransition'+str(args.mg_cutoff)+'\n')
    subjects=open(args.subject_file,'r').read().strip().split('\n')
    num_subjects=str(len(subjects)) 
    transition_dict=dict()
    mg_cutoff=args.mg_cutoff
    start_index=0
    end_index=num_subjects
    if args.start_index!=None:
        start_index=args.start_index
    if args.end_index!=None:
        end_index=min([args.end_index,end_index]) 
    for subject in subjects[start_index:end_index]:
        data=open(subject,'r').read().strip().split('\n')
        compressed_data=compress(data,args.mg_cutoff,args.consecutive_val_thresh)
        transition_dict[subject]=get_number_transition_states(compressed_data)
    for subject in transition_dict:
        #get just the subject id
        filtered_subject=subject.split('/')[-1].split('_')[0] 
        outf.write(filtered_subject+'\t'+str(transition_dict[subject])+'\n')
    print("Done!") 
if __name__=="__main__":
    main()
