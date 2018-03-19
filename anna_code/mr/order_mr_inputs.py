#subjects=open('subjects.txt','r').read().strip().split('\n')
subjects=open('euro.txt','r').read().strip().split('\n')
print("read in subjects")
import sys
source=open(sys.argv[1],'r').read().strip().split('\n') 
outf=open(sys.argv[2],'w')
outf.write(source[0]+'\n')
data_dict=dict()
num_entries=len(source[0].split('\t'))-1
for line in source[1::]:
    tokens=line.split('\t')
    subject=tokens[0]
    data_dict[subject]=line
for subject in subjects:
    if subject in data_dict:
        outf.write(data_dict[subject]+'\n')
    else:
        outf.write(subject+'\t'+'\t'.join(['NA' for i in range(num_entries+1)])+'\n')
        
    
