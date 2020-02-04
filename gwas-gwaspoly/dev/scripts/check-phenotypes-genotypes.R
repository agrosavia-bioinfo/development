#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE)
args = c("agrosavia-phenotype-gwaspoly.tbl", "K5_estructura_tetraploides_2017.txt")

phenotypeFile = args [1]
genotypeFile  = args [2]

phenotype    = read.table (phenotypeFile, header=T, sep=",")
genotype     = read.table (genotypeFile, header=T)

# Matching phenotypes and genotypes
namesPheno = phenotype$snp
namesGeno  = genotype$sample
idxGenos = match (namesGeno, namesPheno)

X = genotype[idxGenos, ]

stopifnot(all(phenotype$snp == rownames(X)))
