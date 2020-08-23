#!/usr/bin/Rscript

source ("lglib.R")

# Add map information (chr, pos) to k-matrix genotype
# INPUTS: a k-matrix genotype, and, map file

#args = c ("agrosavia-genotype-ACGT-CLEANED.tbl",
#		  "agrosavia-genotype-MAP-altref.tbl")

genoFile = args [1]
mapFile  = args [2]

# Read data
geno = read.csv (genoFile, check.names=F)
hd (geno)
map  = read.csv (mapFile)
rownames (map) = map [,1]
hd (map)

# Unify
SNPsNamesGeno = geno [,1]
mapSelected = map [map[,1] %in% SNPsNamesGeno,]
hd (mapSelected)

genoMap = cbind (mapSelected[,c(1,4,5)], geno[,-1])
hd (genoMap)
outFile = addLabel (genoFile, "MAP")

write.table (file=outFile, genoMap, sep=",", quote=F, row.names=F)
		
