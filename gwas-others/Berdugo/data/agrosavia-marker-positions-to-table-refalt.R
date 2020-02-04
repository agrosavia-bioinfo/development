#!/usr/bin/Rscript
# Convert SNPs ref/alt from Berdugo to table

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) 
{
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}
#----------------------------------------------------------
# Util to print head of data
# Fast head, for debug only
#----------------------------------------------------------
hd <- function (data, m=10,n=10) 
{
	msg (deparse (substitute (data)),":")
	if (is.null (dim (data)))
		print (data [1:10])
	else if (ncol (data) < 10) 
		print (data[1:m,])
	else if (nrow (data) < 10)
		print (data[,1:n])
	else 
		print (data [1:m, 1:n])
}

splitSNPs <- function (SNP) {
	refalt = strsplit (SNP, split="/")
	ref    = refalt [[1]][1]
	alt    = refalt [[1]][2]
	return (list (Ref=ref, Alt=alt))
}

data = read.table (file="Positions_markers_chip_8K_papa_jberdugo.txt", header=T)
print (dim (data))
nameMarkers = gsub ("solcap_snp_", "", data$Marker)
refalt = sapply (as.vector (data$SNP), splitSNPs , simplify=T)
colnames (refalt) = nameMarkers
refaltTable = t (refalt)
hd (refaltTable)

message ("Newdata")
newData = cbind (Marker=nameMarkers, refaltTable, Chrom=data [,4], Position=data[,3])
print (dim (newData))
hd (newData)

write.table (file="agrosavia-genotype-map-altref.tbl", newData, row.names=F, quote=F, sep="\t")


