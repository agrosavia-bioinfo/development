#!/usr/bin/Rscript

source ("lglib.R")



#args = c ("agrosavia-genotype-ACGT-CLEANED.tbl",
#		  "agrosavia-genotype-MAP-altref.tbl")

genoFile = args [1]
mapFile  = args [2]

# Read data
geno = read.table (file=genoFile, check.names=F, header=T, sep=",")
hd (geno)
map  = read.table (mapFile, sep="\t", header=T)
rownames (map) = map [,1]
hd (map)

# Unify
SNPsNamesGeno = geno [,1]
mapSelected = map [map[,1] %in% SNPsNamesGeno,]

genoMap = cbind (mapSelected[,c(1,4,5)], geno[,-1])
hd (genoMap)
outFile = addLabel (genoFile, "MAP")

write.table (file=outFile, genoMap, sep=",", quote=F, row.names=F)
		
