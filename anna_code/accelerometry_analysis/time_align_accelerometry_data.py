#align accelerometry data stream by start time and/or day of week
import argparse
from dateutil.parser import parse
from datetime import datetime,timedelta
from Params import *

#Monday = 0, Sunday =6
def get_starting_day(header):
    date_string=header.split('-')[1].strip(' ')
    date=parse(date_string)
    day=date.weekday()
    return day

def parse_args():
    parser=argparse.ArgumentParser(description="align accelerometry data by start time and/or day")
    parser.add_argument("--input_file_list")
    parser.add_argument("--remove_incomplete",default=True) 
    parser.add_argument("--out_prefix")
    return parser.parse_args()
def main():
    args=parse_args()
    input_files=open(args.input_file_list,'r').read().strip().split('\n')
    for fname in input_files:
        #print(str(fname))
        data=open(fname,'r').read().strip().split("\n")
        if args.remove_incomplete==True:
            #check for 120961 datapoints
            if len(data)<120961:
                continue
        #get the starting day.
        starting_day=get_starting_day(data[0])
        #print(str(starting_day))
        observed_order=range(starting_day,7)+range(starting_day)
        #print(str(observed_order))
        day_to_vals=dict()
        data=data[1::]#remove the header 
        for i in range(7):
            key=observed_order[i]
            values=data[start_day_indices[i]:end_day_indices[i]]
            values=[entry.split(',')[0] for entry in values]
            day_to_vals[key]='\n'.join(values)
        output_filename=args.out_prefix+'/'.join(fname.split('/')[-2::])
        #print("output:"+str(output_filename))
        outf=open(output_filename,'w')
        for key in range(7):
            outf.write(day_to_vals[key]+'\n')
        
        
        
if __name__=="__main__":
    main()
    
