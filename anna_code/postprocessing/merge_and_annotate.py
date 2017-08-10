#merges top snp's for each phenotype into a matrix of snp -> phenotype p-values/ phenotype effect sizes
#annotates with maf info & GWAS catalogue
import argparse
from os import listdir
from os.path import isfile,join
import pdb 
def parse_args():
    parser=argparse.ArgumentParser(description="merge and annotate top snp hits for individual traits")
    parser.add_argument("--exclude_list",default=None)
    parser.add_argument("--maf_prefix",default=None)
    parser.add_argument("--gwas_catalog",default=None)
    parser.add_argument("--outf")
    parser.add_argument("--dir_to_merge",help="directory containing the .tsv files to be merged-- all .tsv files in specified directory will be merged") 
    return parser.parse_args()

def get_exclude_dict(fname):
    exclude_dict=dict()
    data=open(fname,'r').read().strip().split('\n')
    for line in data:
        exclude_dict[line]=1
    del data
    return exclude_dict

def get_maf_dict(fname,exclude_dict):
    maf_dict=dict()
    for chrom in range(1,23):
        print("maf for chrom:"+str(chrom))
        data=open(fname+str(chrom)+"_v2.txt",'r').read().split('\n')
        for line in data:
            tokens=line.split()
            if len(tokens)<1:
                print(str(tokens))
            else:
                snp=tokens[0]
                if snp not in exclude_dict: 
                    maf=tokens[-2]
                    maf_dict[snp]=maf
    del data
    return maf_dict 

def get_gc_dict(fname,exclude_snps):
    data=open(fname,'r').read().strip().split('\n')
    gc_dict=dict()
    for line in data[1::]:
        tokens=line.split('\t')
        snp=tokens[2]
        if snp not in exclude_snps: 
            pmid=tokens[-3]
            trait=tokens[-2]
            gc_dict[snp]=[pmid,trait]
    del data
    return gc_dict

def add_task_info(fname,snp_dict):
    data=open(fname,'r').read().strip().split('\n')
    header=data[0].split()
    chrom_index=header.index('chr')
    snp_index=header.index('snp')
    pos_index=header.index('pos')
    a1_index=header.index('a1')
    model_index=header.index('model')
    nmiss_index=header.index('nmiss')
    effect_index=header.index('effect')
    stat_index=header.index('stat')
    pvalue_index=header.index('pvalue')
    trait_index=header.index('trait')
    score_adjustment_index=header.index('score_adjustment')
    ancestry_index=header.index('ancestry')
    region_start_index=header.index('start')
    region_end_index=header.index('end')
    region_size=header.index('size')
    nsnp_index=header.index('n.snps')
    for line in data[1::]:
        tokens=line.split()
        snp=tokens[snp_index]
        #if we are seeing the snp for the first time, add the annotation information 
        if snp not in snp_dict:
            snp_dict[snp]=dict()
            snp_dict[snp]['annotation']=[tokens[chrom_index],snp,tokens[pos_index],tokens[a1_index],tokens[score_adjustment_index],tokens[ancestry_index]]
        #add the trait-specific valuesu for the snp
        snp_dict[snp][tokens[trait_index]]=[tokens[effect_index],tokens[stat_index],tokens[pvalue_index]]
    return snp_dict 
                
def main():
    args=parse_args()
    print("loading annotation files") 
    exclude_dict=get_exclude_dict(args.exclude_list)
    print("loaded dictionary of snps to exclude") 
    maf_dict=get_maf_dict(args.maf_prefix,exclude_dict)
    print("loaded minor allele frequencies")
    gc_dict=get_gc_dict(args.gwas_catalog,exclude_dict)
    print("loaded GWAS catalog annotations") 
    snp_dict=dict()
    tasks=[f for f in listdir(args.dir_to_merge) if isfile(join(args.dir_to_merge,f))]
    for task in tasks:
        print("annotating task:"+str(task))
        full_task_fname=args.dir_to_merge+'/'+task 
        snp_dict=add_task_info(full_task_fname,snp_dict)
        
    #write output!
    print("writing output") 
    outf=open(args.outf,'w')
    tasks=[i.split('.')[0] for i in tasks] 
    outf.write('chr\tsnp\tpos\ta1\tscore_adjustment\tancestry\tmaf\tGWAScatalogPMID\tGWAScatalogTrait\t'+'\t'.join([i+'.effect'+'\t'+i+'.stat'+'\t'+i+'.pvalue' for i in tasks])+'\n')
    for snp in snp_dict:
        if snp in exclude_dict:
            continue #we must exclude this snp, as it's not part of the HRC correctly imputed panel
        outf.write('\t'.join(snp_dict[snp]['annotation']))
        if snp in maf_dict:
            outf.write('\t'+maf_dict[snp])
        else:
            outf.write('\tNA')
        if snp in gc_dict:
            outf.write('\t'+'\t'.join(gc_dict[snp]))
        else:
            outf.write('\t\t') 
        for task in tasks:
            if task in snp_dict[snp]:
                outf.write('\t'+'\t'.join(snp_dict[snp][task]))
            else:
                outf.write('\t'+'\t'*2)
        outf.write('\n')
        
    

    
    

if __name__=="__main__":
    main()
    
    
