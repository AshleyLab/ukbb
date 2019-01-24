#!/bin/bash 
#source: http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=6
#1471 -- atrial fibrillation 


#1081 -- Stroke -  
#1491 -- brain hemmorhage -
#1086 -- subarachnoid hemmorhage 

python get_self_reported_diagnoses.py --diagnosis_field 1471 1081 1086 1491 --outf atrial.fibrillation.1471.tsv 
