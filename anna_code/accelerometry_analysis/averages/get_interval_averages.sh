#!/bin/bash
python get_interval_averages.py --source_files /oak/stanford/groups/euan/projects/ukbb/data/accelerometry/aligned/parts/x$1 \
       --start_interval 21 29 33 41 48 12 18 24 30 \
       --end_interval 29 45 41 48 57 18 24 30 36\
       --outf subject_acceleration_interval_averages.$1.txt
