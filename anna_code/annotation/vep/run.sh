#!/bin/bash
module load perl
module load tabix 
#/home/annashch/ensembl-vep/vep -i vep_input.txt  --port 3337 --plugin PolyPhen_SIFT --cache --force_overwrite
/home/annashch/ensembl-vep/vep -i fitness_vep_input.txt  --port 3337 --plugin PolyPhen_SIFT --cache --force_overwrite 
