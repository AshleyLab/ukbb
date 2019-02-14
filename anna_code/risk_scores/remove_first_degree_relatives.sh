#!/bin/bash
plink2 --bfile afib.cohort.train --remove throw_away.tsv --make-bed --out afib.cohort.train.nofirstdegreerel
plink2 --bfile afib.cohort.test --remove throw_away.tsv  --make-bed --out afib.cohort.test.nofirstdegreerel


