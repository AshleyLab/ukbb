#get the Betas, SE, p-values for exposures and outcomes in physical activity analysis 
import argparse 
from os import listdir
from os.path import isfile, join



def parse_args(): 
    parser=argparse.ArgumentParser(description="get the Betas, SE, p-values for exposures and outcomes in physical activity analysis")
    parser.add_argument("--exposure_files",nargs="+",help="file contains list of variants associated with exposure")
    parser.add_argument("--outcomes",nargs="+",help="string of outcome names")
    parser.add_argument("--exposure_prefix") 
    parser.add_argument("--outcome_prefix") 
    parser.add_argument("--outf") 
    return parser.parse_args() 

def main(): 
    args=parse_args() 
    snp_exposure_dict=dict() 
    for fname in args.exposure_files:
        print(fname) 
        snps=open(fname,'r').read().strip().split('\n')[0].split() 
        for s in snps: 
            if s not in snp_exposure_dict: 
                snp_exposure_dict[s]=dict() 
        full_prefix='/'.join([args.exposure_prefix,fname.split('/')[-1]])
        plink_exposures= [full_prefix+'/'+f for f in listdir(full_prefix) if isfile(join(full_prefix, f))]
        for f_expo in plink_exposures: 
            if ((f_expo.endswith("linear"))or(f_expo.endswith("logistic"))): 
                print('\t'+f_expo)
                data=open(f_expo,'r')
                for line in data: 
                    tokens=line.split() 
                    if tokens[2] in snp_exposure_dict: 
                        if tokens[5] =="ADD": 
                            snp_exposure_dict[tokens[2]][fname.split('/')[-1]]=line 
    print("parsed exposures") 
    #now parse the outcomes 
    for fname in args.outcomes:
        full_prefix='/'.join([args.outcome_prefix,fname])
        print(fname)
        plink_outcomes= [full_prefix+'/'+f for f in listdir(full_prefix) if isfile(join(full_prefix, f))]
        for f_outcome in plink_outcomes: 
            if ((f_outcome.endswith("linear"))or(f_outcome.endswith("logistic"))): 
                print('\t'+f_outcome)
                data=open(f_outcome,'r')
                for line in data: 
                    tokens=line.split() 
                    if tokens[2] in snp_exposure_dict: 
                        if tokens[5] =="ADD": 
                            snp_exposure_dict[tokens[2]][fname]=line 
    print("Writing outputs!")
    outf=open(args.outf,'w')
    for snp in snp_exposure_dict: 
        for fname in snp_exposure_dict[snp]: 
            outf.write(snp+'\t'+fname+'\t'+snp_exposure_dict[snp][fname]+'\n')

if __name__=="__main__": 
    main()
