#!/usr/bin/Rscript

library (dplyr)
source ("lgRlib.R")
#-----------------------------------------------------------
# Select best N SNPs from multiple action models (for GWASpoly and TASSEL)
# PLINK also can produce info of more action models using options
#-----------------------------------------------------------
selectBestN <- function (filename, N) {
	# Read
	data = read.table (file=filename, header=T, sep="\t"); 
	#hd (data,n=13)

	# Select main columns
	dr = data [,c("Marker","GC","MODEL","SCORE", "THRESHOLD", "DIFF")]; 
	#hd (dr,n=13)
	#write.table (file="xdr.tbl", dr, quote=F, sep="\t", row.names=F)

	# Order by N, DIFF, GC
	do = dr [order (dr$MODEL,-dr$DIFF),]; 
	#hd (do,n=15,m=20)
	#write.table (file="xdo.tbl", do, quote=F, sep="\t", row.names=F)

	# Reduce to groups of N
	dm = Reduce (rbind, by(do, do["MODEL"], head, n=N)); 
	#hd (dm, n=13)
	#write.table (file="xdm.tbl", dm, quote=F, sep="\t", row.names=F)

	# Add Count of SNPs
	summ = data.frame (add_count (dm, Marker, sort=T)); 
	hd (summ,n=15,m=20)
	write.table (file="xsumm.tbl", summ, quote=F, sep="\t", row.names=F)

	# Order summary by 
	#best = summ [order (-summ$n, -summ$DIFF),]; 
	best = summ [order (-summ$n, abs (summ$GC - 1)),]; 
	hd (best,n=15,m=20)
	write.table (file="xbest.tbl", best, quote=F, sep="\t", row.names=F)

	# Remove duplicates
	nodups = best [!duplicated (best$Marker, fromLast=F),] ; 
	hd (nodups,n=15,m=20)
	write.table (file="xndups.tbl", nodups, quote=F, row.names=F, sep="\t")

	bestN     = nodups [1:N, ]
	bestModel = nodups [1,"MODEL" ]
	hd (bestN)
	message ("Best model is: ", bestModel)

	# Select SNPs for model and sort by DIFF
	dataModel = data [data[,"MODEL"] %in% bestModel,]
	dataModel = dataModel [order (-dataModel$DIFF),]

	return (dataModel)
}

#filename="tool-TASSEL-scores-Full.csv"
filename="tool-GWASpoly-scores-Full.csv"
bn = selectBestN (filename,10)
hd (bn,m=13)
