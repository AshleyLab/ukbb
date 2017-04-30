#calculates the number of subjects that have data for the fields in the specified *tab file 
import argparse
def parse_args():
    parser=argparse.ArgumentParser("calculates the number of subjects that have data for the fields in the specified *tab file ")
    parser.add_argument("--ukbb_fields_files",nargs="+")
    parser.add_argument("--field_reference_file",default="/scratch/PI/euan/projects/ukbb/file_handlers/ukb_field.tsv")
    parser.add_argument("--outf")
    parser.add_argument("--mode",default="overwrite",help="\"overwrite\" or \"append\" -- the former will overwrite output file specified in --outf; the latter will append to this file")
    return parser.parse_args() 

def main():
    args=parse_args()
    if args.mode=="overwrite":
        outf=open(args.outf,'w')
    elif args.mode=="append":
        outf=open(args.outf,'a')
    else:
        exit("invalid --mode argument:"+args.mode+" was specified; must be one of \"overwrite\" or \"append\"")
    #build a dictionary of field id to field & category name
    ukb_fields=open(args.field_reference_file,'r').read().strip().split('\n')
    name_map=dict()
    for line in ukb_fields:
        tokens=line.split('\t')
        category_id=tokens[0]
        field_id=tokens[1]
        category_name=tokens[2]
        field_name=tokens[3]
        name_map[field_id]=[field_id,category_id,category_name,field_name]
    #print(str(name_map))
    print("built dictionary mapping field id to field name!")
    
    field_dict=dict() 
    for data_file in args.ukbb_fields_files:
        print(str(data_file))
        counter=0
        got_header=False
        with open(data_file,'r') as f:
            for line in f:
                if counter%100000==0:
                    print(str(counter))
                if got_header==False:
                    header=line.split('\t')
                    for field in header[1::]:
                        field_dict[field]=0
                    num_tokens=len(header) 
                    got_header=True
                else:
                    tokens=line.split('\t')
                    for i in range(1,num_tokens):
                        cur_field=header[i]
                        cur_val=tokens[i]
                        if cur_val!="NA":
                            field_dict[cur_field]+=1
                counter+=1
    print("built field id to field name mapping") 
    #write the output file
    print("generating output file!") 
    for field in field_dict:
        #get the field name & category information
        try:
            metadata=name_map[field.split('.')[1]]
        except:
            print('WARNING!:field'+str(field)+' is not found in the field name set -- you will need to look it up manually')
            metadata=[field,'','','']
        outf.write('\t'.join(metadata)+'\t'+str(field_dict[field])+'\n')
        
    
if __name__=="__main__":
    main()
    
