#!/bin/bash

GENOTYPE=genob

### Step 1 ### 

# Investigate missingness per individual and per SNP and make histograms.
plink --bfile $GENOTYPE --missing    
# output: plink.imiss and plink.lmiss, these files show respectively the proportion of missing SNPs per individual and the proportion of missing individuals per SNP.

# Generate plots to visualize the missingness results.
Rscript --no-save scripts/hist_miss.R

