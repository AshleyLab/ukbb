#assemble the phenotype file for a GWAS analysis of physical activity
#Need to avoid reading in full files to prevent memory errors 
import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="assemble fields for a GWAS analysis of physical activity")
    parser.add_argument("--tables",nargs='+')
    parser.add_argument("--fields",nargs='+')
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    subject_dict_continuous=dict()
    subject_dict_categorical=dict() 
    traits_continuous=[]
    traits_categorical=[] 
    #get the bulk data tables loaded first
    for i in range(len(args.tables)):
        table=open(args.tables[i],'r')
        header=table.readline().strip('\n').split('\t')
        entries=open(args.fields[i],'r').read().strip().split('\n')
        fields=[f.split('\t')[0] for f in entries]
        field_names=[f.split('\t')[1] for f in entries]
        field_type=[f.split('\t')[2] for f in entries]
        for j in range(len(field_names)):
            if field_type[j]=="categorical":
                traits_categorical.append(field_names[j])
            else:
                traits_continuous.append(field_names[j])
        print("read in UKBB table:"+args.tables[i])
        #get the index of each field
        field_indices=[]
        for field in fields:
            field_indices.append(header.index(field))
        line=header
        while line!="":
            line=table.readline().strip('\n')
            if line=="":
                break
            tokens=line.split('\t')
            subject=tokens[0]
            if subject not in subject_dict_continuous:
                subject_dict_continuous[subject]=dict()
                subject_dict_categorical[subject]=dict() 
            for j in range(len(field_indices)):
                value=tokens[field_indices[j]]
                if value=="NA":
                    value="-1000" 
                if field_type[j]=="categorical": 
                    subject_dict_categorical[subject][field_names[j]]=value
                else:
                    subject_dict_continuous[subject][field_names[j]]=value                
    #generate an output file of phenotypes
    subjects=subject_dict_continuous.keys()
    outf_continuous=open(args.outf+".continuous.txt",'w')
    outf_categorical=open(args.outf+".categorical.txt",'w')
    outf_continuous.write('FID\tIID\t'+'\t'.join(traits_continuous)+'\n')
    outf_categorical.write('FID\tIID\t'+'\t'.join(traits_categorical)+'\n')
    for subject in subjects:
        out_string_continuous=subject+'\t'+subject
        out_string_categorical=subject+'\t'+subject
        for trait in traits_continuous:
            if trait in subject_dict_continuous[subject]:
                out_string_continuous=out_string_continuous+'\t'+subject_dict_continuous[subject][trait]
            else:
                out_string_continuous=out_string_continuous+'\t'+'-1000'
        outf_continuous.write(out_string_continuous+'\n')
        for trait in traits_categorical:
            if trait in subject_dict_categorical[subject]:
                out_string_categorical=out_string_categorical+'\t'+subject_dict_categorical[subject][trait]
            else:
                out_string_categorical=out_string_categorical+'\t'+'-1000'
        outf_categorical.write(out_string_categorical+'\n')

if __name__=="__main__":
    main() 
    
