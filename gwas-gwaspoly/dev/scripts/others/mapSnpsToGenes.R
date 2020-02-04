#!/usr/bin/Rscript

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) {
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}

args = commandArgs(trailingOnly = TRUE)
options (width=300)

args = c ("snps.tbl", "genes.gff")
snpsFile  = args [1]
genesFile = args [2]


if (file.exists ("map.RData")) {
	load (file="map.RData")
}else{ 
	msg ("Reading snps...")
	#snpsTable  = read.table (file=snpsFile, header=T)
	msg ("Reading genes...")
	genesTable = read.table (file=genesFile, header=F, sep="\t")
}

snpsTable  = read.table (file=snpsFile, header=T)
print (snpsTable)
n = length (snpsTable)

mapTableMat = matrix (ncol=5)
#colnames (mapTableMat) = c("snp", "chr", "pos", "type", "id")
mapTable = as.table (mapTableMat)

for (i in 1:n) {
	snp    = snpsTable [i,]
	idsnp  = snp [4]
	chr    = sprintf ("ST4.03ch%02d",strtoi (snp[5]))
	pos    = snp [6]
	scr    = snp [9]
	msg (chr, pos, scr)


	genes = genesTable [genesTable$V4>=pos & genesTable$V4<=pos & genesTable$V1==chr,]
	print (dim(genes))
	if (nrow(genes) > 0){
		prd = genes [,c(3,9)]
		map = cbind (c(idsnp), c(chr), c(pos), prd)
		mapTableMat = rbind (mapTableMat, map)
		#mapTable = rbind (mapTable, genes[,c(9)])
		#print (mapTable)
		quit()
	}
}



#save(snpsTable, genesTable, file="map.RData") 



