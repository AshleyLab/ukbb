import argparse 
def parse_args(): 
    parser=argparse.ArgumentParser(description="collapses a matrix of ICD codes to the specified number of characters")
    parser.add_argument("--input_icd_file") 
    parser.add_argument("--num_chars",type=int,default=3) 
    parser.add_argument("--outf") 
    return parser.parse_args() 

def main(): 
    args=parse_args() 

    #map patient->collapsed_icd->1 
    patient_to_collapsed_icd=dict()
    
    #set of collapsed icd's
    collapsed_icd=set([])

    data=open(args.input_icd_file,'r')

    #map column index -> collapsed icd it corresponds with 
    col_to_collapsed_icd=dict()
    line=data.readline()
    header=line.split('\t')  

    for i in range(1,len(header)): 
        cur_collapsed_icd=header[i][0:args.num_chars]
        col_to_collapsed_icd[i]=cur_collapsed_icd 
        collapsed_icd.add(cur_collapsed_icd) 
    counter=1
    while line!="": 
        line=data.readline() 
        tokens=line.split('\t') 
        subject=tokens[0] 
        patient_to_collapsed_icd[subject]=dict() 
        for i in range(1,len(tokens)): 
            if tokens[i]=="1": 
                cur_icd=col_to_collapsed_icd[i]
                patient_to_collapsed_icd[subject][cur_icd]="1"
        line=data.readline() 
        counter+=1 
        if counter%1000==0: 
            print(str(counter))

    print("writing output file") 
    outf=open(args.outf,'w') 
    collapsed_icd=list(collapsed_icd) 
    outf.write("IID\t"+"\t".join(collapsed_icd)+'\n')
    for patient in patient_to_collapsed_icd: 
        outf.write(patient) 
        for icd in collapsed_icd: 
            if icd in patient_to_collapsed_icd[patient]: 
                outf.write('\t1') 
            else: 
                outf.write('\t0') 
        outf.write('\n') 

if __name__=="__main__": 
    main() 
