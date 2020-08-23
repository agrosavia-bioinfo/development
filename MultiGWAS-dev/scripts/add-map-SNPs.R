#!/usr/bin/Rscript

source ("lglib.R")

# Add map information (chr, pos) to k-matrix genotype

#args = c ("agrosavia-genotype-ACGT-CLEANED.tbl",
#		  "agrosavia-genotype-MAP-altref.tbl")

SNPsFile = args [1]
mapFile  = args [2]

# Read data
SNPs = read.csv (SNPsFile, check.names=F)
hd (SNPs)
map  = read.csv (mapFile)
rownames (map) = map [,1]
hd (map)

# Unify
SNPsNames = SNPs [,1]
mapSelected = map [map [,1] %in% SNPsNames,]
mapSelected = mapSelected [match (SNPsNames, mapSelected [,1]),]
hd (mapSelected)

SNPsMap = cbind (SNPs [, c(1,2,3)], mapSelected[,c(4,5)])
hd (SNPsMap)
outFile = addLabel (SNPsFile, "MAP")

write.table (file=outFile, SNPsMap, sep=",", quote=F, row.names=F)
		
