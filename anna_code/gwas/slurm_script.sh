#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
export chrom=chrom$1
#################
#set a job name
SBATCH --job-name=PLINK_$chrom
#################
#a file for job output, you can check job progress, append the job ID with %j to make it unique
SBATCH --output=PLINK_$chrom.out 
#################
# a file for errors from the job
SBATCH --error=PLINK_$chrom.err 
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
SBATCH --time=5:00
#################
#Quality of Service (QOS); think of it as job priority, there is also --qos=long for with a max job length of 7 days, qos normal is 48 hours.
# REMOVE "normal" and set to "long" if you want your job to run longer than 48 hours,
# NOTE- in the hns partition the default max run time is 7 days , so you wont need to include qos

#SBATCH --qos=normal

# We are submitting to the dev partition, there are several on sherlock: normal, gpu, owners, hns, bigmem (jobs requiring >64Gigs RAM)
#
SBATCH -p euan,owners
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#memory per node; default is 4000 MB per CPU
#SBATCH --mem=4000
#################

# Have SLURM send you an email when the job ends or fails, careful, the email could end up in your clutter folder
# Also, if you submit hundreds of jobs at once you will get hundreds of emails.

#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=YourSUNetID@stanford.edu


#now run normal batch commands
# note the "CMD BATCH is an R specific command
#module load R/3.3.0
# You can use srun if your job is parallel
#srun R CMD BATCH  ./rtest.R

# otherwise:

#R CMD BATCH  ./rtest.R
plink --bfile /scratch/PI/euan/projects/ukbb/data/genetic_data/calls/plink/$chrom --pheno accelerometry.pheno --assoc  --all-pheno --pfilter 1e-6 --missing-phenotype -1000 --out $chrom/$chrom
