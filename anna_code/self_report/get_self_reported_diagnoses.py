import argparse #python library to read in arguments 
def parse_args(): 
    parser=argparse.ArgumentParser() 
    # nargs="+" lets you give a list of values for this field, separated by spaces 
    parser.add_argument("--diagnosis_fields",nargs="+", help="http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=6")
    parser.add_argument("--source_table",default="/oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb7454.tab")
    parser.add_argument("--ukbb_code",default="f.20002")
    parser.add_argument("--outf")
    return parser.parse_args() 

def main(): 
    args=parse_args() 
    #1. create a dictionary of the diagnosis_fields we are interested in (for fast lookups) 
    diagnosis_fields=dict() 
    for field in args.diagnosis_fields: 
        diagnosis_fields[field]=1 
    print("populated dictionary with fields of interest") 
    all_fields=list(diagnosis_fields.keys()) 

    #open the source table in read-only mode, as denoted by 'r'. So we can read from the file, but not modify it 
    data=open(args.source_table,'r').read().strip().split('\n')
    print("read in source table") 
    header=data[0].split('\t')  
    
    #split the header line to get indices of f.20002
    columns=[] 
    for i in range(len(header)): 
        if header[i].startswith(args.ukbb_code): 
            columns.append(i) 
    print("indexed header") 
    
    patient_dict=dict() 
    for line in data[1::]:
        #convert the string of text to a list 
        tokens=line.split('\t') 
        patient=tokens[0] 
        #look at the columns that contain f.20002 entries 
        for c in columns: 
            #check to see if column value is an "allowed" diagnosis 
            if tokens[c] in diagnosis_fields: 
                #if it is, store the patient<-> diagnosis mapping in dictionary 
                if patient not in patient_dict: 
                    patient_dict[patient]=dict() 
                patient_dict[patient][tokens[c]]=1

    print("finishing parsing file") 
    outf=open(args.outf,'w') 
    #writing header 
    outf.write("Patient"+'\t'+'\t'.join(all_fields)+'\n')
    for patient in patient_dict: 
        outf.write(patient) 
        for cur_field in all_fields: 
            if cur_field in patient_dict[patient]: 
                outf.write('\t1') 
            else: 
                outf.write('\t0') 
        outf.write('\n') 

if __name__=="__main__": 
    main() 

