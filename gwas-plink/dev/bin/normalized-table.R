#!/usr/bin/Rscript

normalize <- function(x) { return ((x - min(x)) / (max(x) - min(x))) }

args = commandArgs(trailingOnly = TRUE)
args = c("phenotype-plink.tbl")


filename    = args [1]
outFilename = paste0(strsplit (filename, split="[.]")[[1]][1],"-norm.tbl")

data   = read.table (file=filename, header=T)
data[,3] = normalize (data[,3])

write.table (file = outFilename, data, quote=F, row.names=F)
