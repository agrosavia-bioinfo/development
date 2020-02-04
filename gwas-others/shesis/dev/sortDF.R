#!/usr/bin/Rscript

suppressMessages (library (dplyr))
options (width=300, stringsAsFactors=T)
args = commandArgs(trailingOnly=T)

args = c("out-shesis-adj-AGs.txt", "8")

filename = args [1]
numCol   = as.integer (args [2])

data = read.table (file=filename, header=T )
sdata = data [order (data[,numCol]),]

name = strsplit (filename, split="[.]")[[1]][1]
ext  = strsplit (filename, split="[.]")[[1]][2]

outName = paste0 (name, "-sorted.", ext)

write.table (file=outName, sdata, quote=F, row.names=F, sep="\t")

