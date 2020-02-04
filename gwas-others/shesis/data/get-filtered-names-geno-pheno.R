#!/usr/bin/Rscript

# Get sample and marker names from filtered genotype and phenotype

suppressMessages (library (dplyr))
options (width=300, stringsAsFactors=T)
args = commandArgs(trailingOnly=T)

args = c ("geno-6-FLT-gwasp.tbl","pheno-FLT-gwasp.tbl")
genotype = args [1]
phenotype = args [2]

gNames = read.csv (file=genotype)
pNames = read.csv (file=phenotype)

write.table (file="markers-filtered-names.tbl", gNames[,1], col.names=c("Markers"), quote=F,row.names=F)
write.table (file="samples-filtered-names.tbl", paste0("A", pNames[,1]), col.names=c("Samples"), quote=F,row.names=F)
