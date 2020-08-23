#!/usr/bin/Rscript



hd <- function (data, m=10,n=10) {
	message (deparse (substitute (data)),":")
	if (is.null (dim (data)))
		print (data [1:10])
	else if (ncol (data) < 10) 
		print (data[1:m,])
	else if (nrow (data) < 10)
		print (data[,1:n])
	else 
		print (data [1:m, 1:n])
}
args = commandArgs(trailingOnly = TRUE)
args = c ("plink")

name = args [1]
pedFile = paste0 (name, ".ped")
mapFile = paste0 (name, ".map")

ped = read.table (file=pedFile, header=F, sep="\t")
hd (ped)
map = read.table (file=mapFile, header=F, sep="\t")
hd (map)

namesInd = ped [,2]
hd (namesInd)

snps  = t (ped [, -c(1:6)])
colnames (snps) = namesInd
hd (snps)

geno = data.frame (SNP=map[,2],CHR=map[,1],POS=map[,4], snps)
hd (geno)

genoFilename = paste0 (name, ".tbl")
write.table (file=genoFilename, geno, sep="\t", quote=F, rownames=F)
