import argparse
from os import listdir
from os.path import isfile,join
import numpy as np

def parse_args():
    parser=argparse.ArgumentParser(description= "get average acceleration values for specified intervals of time")
    parser.add_argument("--source_files",help="file that lists all subject source files")
    parser.add_argument("--start_interval",help="start values (in hours on 0 - 24 scale) for intervals",type=int,nargs="+")
    parser.add_argument("--end_interval",help="end values (in hours on 0 - 24 scale) for intervals", type=int, nargs="+")
    parser.add_argument("--outf",help="output file")
    parser.add_argument("--data_start",help="first datapoint, in seconds since 12:00:00am",default=7200,type=float)
    parser.add_argument("--sample_rate",type=float,default=12,help="samples rate per minute") 
    return parser.parse_args()

def format_intervals(args):
    
    '''
    Obtain the first interval start position, 
    the number of points in the interval, 
    and the space between adacent intervals 
    '''
    
    interval_dict=dict()
    start_interval=args.start_interval
    end_interval=args.end_interval
    sample_rate=args.sample_rate
    data_start=args.data_start
    
    for i in range(len(start_interval)):
        first_start=start_interval[i]*60*sample_rate - data_start
        interval_length=(end_interval[i]-start_interval[i])*60*sample_rate
        
        #get the gap to the next interval
        interval_dict[i]=[first_start,interval_length]        
    return interval_dict 

def get_subject_interval_mean(data,data_dict,interval_dict,subject,sample_rate):
    num_points=data.shape[0]
    data_dict[subject]=dict() 
    for i in interval_dict:
        start_val=interval_dict[i][0]
        interval_length=interval_dict[i][1]
        end_val=start_val+interval_length
        cur_interval=data[start_val:end_val,0]
        data_dict[subject][i]=cur_interval
        start_val=start_val+24*60*sample_rate 
        end_val=start_val+interval_length 
        while end_val < num_points:
            cur_interval=data[start_val:end_val,0]
            data_dict[subject][i]=np.concatenate((cur_interval,data_dict[subject][i]))
            start_val=start_val+24*60*sample_rate
            end_val=start_val+interval_length
            #get the average value
        data_dict[subject][i]=np.nanmean(data_dict[subject][i])
    return data_dict


def main():
    #read in the arguments 
    args=parse_args()
    data_dict=dict() # subject -> interval_start -> average value
    interval_dict=format_intervals(args)
    subject_files=open(args.source_files,'r').read().strip().split('\n')
    for f in subject_files: 
        #load in the data
        print(str(f))
        data=np.genfromtxt(f,skip_header=1,delimiter=',')
        subject=f.split('/')[-1].split('_')[0]
        #get the interval averages for the subject 
        data_dict=get_subject_interval_mean(data,data_dict,interval_dict,subject,args.sample_rate)
                
    #generate the output file
    outf=open(args.outf,'w')

    #get the header:
    outf.write('Subject')
    interval_indices=list(interval_dict.keys())
    interval_indices.sort() 
    for i in interval_indices: 
        cur_interval=str(args.start_interval[i])+'_'+str(args.end_interval[i])
        outf.write('\t'+cur_interval)
    outf.write('\n')

    #write subject interval values
    for subject in data_dict:
        outf.write(subject)
        for i in interval_indices:
            outf.write('\t'+str(data_dict[subject][i]))
        outf.write('\n')

        
                        

                        
                







                
                
    

if __name__=="__main__":
    main()
    
