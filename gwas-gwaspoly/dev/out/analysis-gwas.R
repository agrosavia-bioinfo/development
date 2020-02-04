#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE)
args = c("tetra", "diplo")

diploDir = args [1]
tetraDir = args [2]

diploFiles = paste0(diploDir, "/",list.files (path=diploDir, pattern=".tbl"))
print (diploFiles)

f1 = diploFiles [1]
print (f1)
d1 = read.csv (file=f1, header=T, sep="\t")

snpsTable = data.frame ()
snpsTable = rbind (snpsTable, data.frame (d1[,1],d1[,2],d1[,3]))
write.table (file="f1.tbl", snpsTable)


