import argparse
import pdb

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument('--motif_hits')
    parser.add_argument('--vep_annotation')
    parser.add_argument('--outf')
    return parser.parse_args()

def make_vep_dict(vep_file):
    data=open(vep_file,'r').read().strip().split('\n')
    vep_dict=dict()
    for line in data:
        if line.startswith('#'):
            continue
        tokens=line.split('\t')
        location=tokens[1].split(':')
        chrom='chr'+location[0]
        pos=location[1].split('-')[0] 
        var_type=tokens[3]
        impact=tokens[4]
        gene=tokens[5]
        function=tokens[9]
        entry=tuple([chrom,pos])
        vep_dict[entry]=[var_type,impact,gene,function] 
    return vep_dict
    
def main():
    args=parse_args()
    vep_dict=make_vep_dict(args.vep_annotation)
    print("made vep dict")
    motif_hits=open(args.motif_hits,'r').read().strip().split('\n')
    outf=open(args.outf,'w')
    outf.write(motif_hits[0]+'\t'+'Consequence\tImpact\tGene\tBiotype\n')
    for line in motif_hits[1::]:
        tokens=line.split('\t')
        chrom=tokens[0]
        pos=tokens[1]
        entry=tuple([chrom,pos])
        if entry not in vep_dict:
            print(str(entry))
            continue
        vep_info=vep_dict[entry]
        outf.write(line+'\t'+'\t'.join(vep_info)+'\n')
        
if __name__=="__main__":
    main()
    
    
