#!/usr/bin/Rscript

#options (width=300, stringsAsFactors=F)
#args = commandArgs(trailingOnly=T)
library (parallel)
suppressMessages (library (dplyr))

#------------------------------------------------------------------------------
## Format gwaspoly phenotype to plink format
#------------------------------------------------------------------------------
gwaspPhenoToPlinkPheno <- function (gwaspPhenoFile, outFile="") 
{
	message ("Creating plink phenotype...")
	phenotype = read.csv (file=gwaspPhenoFile, header=T)
	idNames = as.character (phenotype [,1])

	samples = phenotype [,1]
	traits  = phenotype [,2]

	message (">>> Writing plink phenotype...")
	plinkPheno = data.frame (FID=0,IID=samples, TRAIT= traits)
	if (outFile=="")
		outFile = paste0 ("out/",strsplit (gwaspPhenoFile, split="[.]")[[1]][1], "-plink.tbl")
	write.table (file=outFile, plinkPheno, col.names=T, row.names=F, quote=F, sep="\t")
}

#------------------------------------------------------------------------------
## Create plink MAP file from gwaspoly genotype 
#------------------------------------------------------------------------------
gwaspGenoToPlinkMap <- function (gwaspGenoFile) 
{
	message (">>> Creating plink MAP file...")
	genotype    <<- read.table (file=gwaspGenoFile, header=T,stringsAsFactors=T,sep=",")
	markers     <<- as.character(genotype [,1])
	chromosomes <<- genotype [,2]
	positions   <<- genotype [,3]

	plinkMap     <<- data.frame (chr=chromosomes, iid=markers, dist=0, positions=positions)
	plinkMapSorted <<- plinkMap %>% arrange (chr, positions)
	outFile   = paste0 ("out/",strsplit (gwaspGenoFile, split="[.]")[[1]][1], "-plink.map")
	write.table (file=outFile, plinkMapSorted, col.names=F, row.names=F, quote=F, sep="\t")
	return (plinkMapSorted$iid)
}


#----------------------------------------------------------
# Getting ref/alt alleles for SNPs
#----------------------------------------------------------
bases <- c("A","C","G","T")
get.ref <- function(x) {
	y <- paste(na.omit(x),collapse="")
	ans <- apply(array(bases),1,function(z,y){length(grep(z,y,fixed=T))},y)
	if (sum(ans)>2) {stop("Error in genotype matrix: More than 2 alleles")}
	if (sum(ans)==2) {ref.alt <- bases[which(ans==1)]}
	if (sum(ans)==1) {ref.alt <- c(bases[which(ans==1)],NA)}

	return(ref.alt)
}
#----------------------------------------------------------
# Create plink PED file from gwaspoly genotype 
#----------------------------------------------------------
gwaspTetraGenoToPlinkPed <- function (gwaspGenoFile, markersIdsMap) 
{
	plinkFile  = paste0 ("out/", strsplit (gwaspGenoFile, split="[.]")[[1]][1], "-plink")

	if (file.exists (paste0(plinkFile,".ped"))) {
		msg ("Loading plink file...")
		return (plinkFile)
	}else {
		message (">>> Creating plink PED file...")
		genotype   = read.csv (file=gwaspGenoFile, header=T)
		alleles    <- as.matrix (genotype [,-c(1,2,3)])
		rownames (alleles) <- genotype [,1]

		message (">>> Creating transposed genotype...")
		markersIds        <- genotype [,1] 
		samplesIds        <- colnames (alleles)

		msg ("Getting Ref/Alt Alleles...")
		refAltAlleles <- apply(alleles,1,get.ref)

		# Convert tetra to diplo 
		allelesDiplo  <- tetraToDiplos (alleles, refAltAlleles)
		rownames (allelesDiplo) = markersIds
		colnames (allelesDiplo) = samplesIds

		# Adjust for plink PED file
		allelesPlink <- t(allelesDiplo[markersIdsMap,])
		genoPED    <- cbind (0, samplesIds, 0,0,0,-9, allelesPlink)

		msg ("Writing plink diplo PED file...")
		plinkFilePed  = paste0 (plinkFile, ".ped")
		write.table (file=plinkFilePed, genoPED, col.names=F, row.names=F, quote=F, sep="\t")

		# Other files
		# Write tetra to diplo for gwasp
		gwasp2genotype = cbind (genotype [,1:3], allelesDiplo)
		outFile   = paste0 ("out/",strsplit (gwaspGenoFile, split="[.]")[[1]][1], "-diplo.tbl")
		write.table (file=outFile, gwasp2genotype, col.names=F, row.names=F, quote=F, sep="\t")

		# Transposed gwasp file
		transposedAlleles <- t(alleles)
		rownames (transposedAlleles) = samplesIds
		colnames (transposedAlleles) = markersIds
		write.table (file ="out/tmp-plink-transposed.tbl", transposedAlleles, col.names=T, row.names=T, quote=F, sep="\t")


		# Write diplo matrix for original matrix
		write.table (file="out/tmp-diplosMatrix.tbl",allelesDiplo,  quote=F, sep="\t")
	}
	
	return (plinkFile)
}
#----------------------------------------------------------
# Add tabs to alleels changign AG --> A	G
#----------------------------------------------------------
tetraToDiplos <- function (alleles, refAltAlleles) {
	#alleles = t (alleles)
	msg ("Converting tetra to diplos...")
	matList <- apply (alleles, 2, list) # Lists by SNPs

	t2d <- function (allele, ref, alt, snp) {
		if (allele=="0000") return ("00")
		else if (strrep (ref,4) == allele) return (paste0(ref,ref))
		else if (strrep (alt,4) == allele) return (paste0(alt,alt))
		else return (paste0(ref,alt))
	}
	ref <- refAltAlleles [1,]
	alt <- refAltAlleles [2,]

	diplosLst <<- list()
	for (ls in matList) {
		diplos    <- mcmapply (t2d,ls[[1]], ref, alt, rownames (alleles), mc.cores=4)
		diplosLst <- append (diplosLst, diplos)
	}

	allelesDiplo <- matrix (unlist(diplosLst), ncol=ncol(alleles), nrow=nrow(alleles), byrow=F)

	return (allelesDiplo)
}

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) 
{
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}
#----------------------------------------------------------
# Main
#----------------------------------------------------------
main <- function () {
	args = c ("agrosavia-genotype-tetra-ACGT.tbl", "agrosavia-phenotype.tbl")

	gwaspGenoFile  = args [1]
	gwaspPhenoFile = args [2]

	message (">>> Converting gwaspoly to plink formats...")
	gwaspPhenoToPlinkPheno (gwaspPhenoFile)
	markersIdsMap = gwaspGenoToPlinkMap (gwaspGenoFile)
	gwaspTetraGenoToPlinkPed (gwaspGenoFile, markersIdsMap)
}

