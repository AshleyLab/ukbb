#identifies snp's below a p-value in discovery and replication cohort.
import argparse
from os import listdir
from os.path import isfile,join

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--discovery_dir")
    parser.add_argument("--replication_dir")
    parser.add_argument("--pval_discovery_thresh",type=float,default=5e-8)
    parser.add_argument("--pval_validation_thresh",type=float,default=0.01) 
    parser.add_argument("--outf")
    parser.add_argument("--features",nargs="+")
    parser.add_argument("--bidirectional",action='store_true',default=False)
    parser.add_argument("--directly_genotyped",action='store_true',default=False,help="consider only directly genotyped SNPs in the analysis")
    parser.add_argument("--directly_genotyped_metadata",default=None)
    parser.add_argument("--maf_thresh",action="store_true",default=False)
    parser.add_argument("--maf_min",type=float,default=None,help="minimum maf to consider")
    parser.add_argument("--maf_max",type=float,default=None,help="maximum maf to consider")
    parser.add_argument("--maf_metadata")
    return parser.parse_args()

def get_call_dict(calls_dir):
    calls=dict()
    for chrom in range(1,23):
        print("parsing calls for chrom:"+str(chrom))
        call_data=open(calls_dir+"/"+"ukb_snp_chr"+str(chrom)+"_v2.bim",'r').read().strip().split('\n')
        for line in call_data:
            tokens=line.split('\t')
            chrom=tokens[0]
            pos=tokens[3]
            entry=tuple([chrom,pos])
            calls[entry]=1
    print("generated dictionary of directly genotyped SNPs")
    print(str(len(calls.keys())))
    return calls
def get_maf_dict(maf_dir,maf_min,maf_max):
    maf_dict=dict()
    for chrom in range(1,23):
        print("parsing maf for chrom:"+str(chrom))
        maf_data=open(maf_dir+"/"+"ukb_mfi_chr"+str(chrom)+"_v2.txt",'r').read().strip().split('\n')
        for line in maf_data:
            tokens=line.split()
            if len(tokens)<1:
                print(str(tokens))
            else:
                snp=tokens[0]
                maf=float(tokens[-2])
                if maf > maf_min:
                    if maf < maf_max:
                        maf_dict[snp]=maf
    print("generated minor allele dictionary")
    return maf_dict


def main():
    args=parse_args()
    if (args.directly_genotyped==True):
        call_dict=get_call_dict(args.directly_genotyped_metadata)
    if (args.maf_thresh==True):
        maf_dict=get_maf_dict(args.maf_metadata,args.maf_min,args.maf_max) 
    discovery_features=args.features 
    outf_summary=open(args.outf+"."+"summary",'w')
    outf_discovery_bed=open(args.outf+"."+"discovery.bed",'w')
    outf_replication_bed=open(args.outf+"."+"replication.bed",'w')
    for feature in discovery_features:
        outf=open(args.outf+"."+feature,'w')
        for chrom in range(1,23):
            discovery_dict=dict()
            try:
                cur_file=open(args.discovery_dir+'/'+feature+'/'+feature+'.'+str(chrom)+'.'+'continuous.assoc.linear').read().strip().split('\n')
            except:
                cur_file=open(args.discovery_dir+'/'+feature+'/'+feature+'.'+str(chrom)+'.'+'categorical.qassoc').read().strip().split('\n')
            for line in cur_file[1::]:
                tokens=line.split()
                if len(tokens)<3:
                    print(str(tokens))
                    continue 
                pval=tokens[-1]
                try:
                    pval=float(pval)
                except:
                    continue 
                rs=tokens[1]
                chrom=tokens[0]
                pos=tokens[2]
                entry=tuple([chrom,pos])
                if ((args.directly_genotyped==True) and (entry not in call_dict)):
                    continue
                if ((args.maf_thresh==True) and (rs not in maf_dict)):
                    continue 
                if pval<=args.pval_discovery_thresh:
                    discovery_dict[rs]=tokens
            print("finished discovery set for chrom:"+str(chrom))
            replication_dict=dict()
            if args.replication_dir!=None:
                try:
                    cur_file=open(args.replication_dir+'/'+feature+'/'+feature+'.'+str(chrom)+'.'+'continuous.assoc.linear').read().strip().split('\n')
                except:
                    cur_file=open(args.replication_dir+'/'+feature+'/'+feature+'.'+str(chrom)+'.'+'categorical.qassoc').read().strip().split('\n')
                for line in cur_file[1::]:
                    tokens=line.split()
                    if len(tokens)<3:
                        print(str(tokens))
                        continue 
                    pval=tokens[-1]
                    try:
                        pval=float(pval)
                    except:
                        continue 
                    rs=tokens[1]
                    chrom=tokens[0]
                    pos=tokens[2]
                    entry=tuple([chrom,pos])
                    if ((args.directly_genotyped==True) and (entry not in call_dict)):
                        continue
                    if ((args.maf_thresh==True) and (rs not in maf_dict)):
                        continue 
                    if pval<=args.pval_validation_thresh:
                        replication_dict[rs]=tokens
                print("finished replication set for chrom:"+str(chrom))            
            num_discovered=0
            num_replicated=0
            if args.bidirectional==True: 
                num_new=0 
            for snp in discovery_dict:
                num_discovered+=1
                outf.write('\t'.join(discovery_dict[snp]))
                if snp in replication_dict:
                    num_replicated+=1
                    outf.write('\t'+'\t'.join(replication_dict[snp]))
                outf.write('\n')
            if args.bidirectional==True:
                if len(discovery_dict.values())>0: 
                    num_to_skip=len(discovery_dict.values()[0])
                else:
                    num_to_skip=9 
                for snp in replication_dict:
                    if snp not in discovery_dict:
                        snp_pval=float(replication_dict[snp][-1])
                        if snp_pval <= args.pval_discovery_thresh:
                            outf.write('\t'*num_to_skip+'\t'.join(replication_dict[snp])+'\n')
                            num_new+=1 
            outf_summary.write(feature+'\t'+str(chrom)+'\t'+str(num_discovered)+'\t'+str(num_replicated))
            if args.bidirectional==True:
                outf_summary.write('\t'+str(num_new))
            outf_summary.write('\n')
            
if __name__=="__main__":
    main()
    
