original_data=open("timed_out.redo").read().strip().split('\n')
#redo_data=open("timed_out.redo").read().strip().split('\n')
new_data=open("timed_out.redo2").read().strip().split('\n')

#new-data - original data + redo_data
original_dict=dict()
for line in original_data:
    original_dict[line]=1

#redo_dict=dict()
#for line in redo_data:
#    redo_dict[line]=1

new_dict=dict()
for line in new_data:
    new_dict[line]=1

keep=[]
for line in new_data:
    if line in original_dict:
        continue
    else:
        keep.append(line)

#keep=keep+redo_data
outf=open('net_set.txt','w')
outf.write('\n'.join(keep))


