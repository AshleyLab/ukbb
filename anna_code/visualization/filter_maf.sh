#!/bin/bash
python filter_maf.py --outf $1.filtered \
       --maf_filter 0.01 \
       --maf_pickle_prefix /oak/stanford/groups/euan/projects/ukbb/code/anna_code/useful_pickles/maf_ \
       --inputf $1.aggregate.maf.0.01
