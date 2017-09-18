#!/bin/sh
$HOME/apps/plink/plink --bfile /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/plink_dir_genotype/ukb_cal_chr$2\_v2 --clump overall_acceleration_average_inputs_not_normalized/$1/$1.$2.continuous.$1.glm.linear --out overall_acceleration_average_inputs_not_normalized/$1/$1.$2
