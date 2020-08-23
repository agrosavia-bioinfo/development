#!/usr/bin/Rscript

#-------------------------------------------------------------
# Return the format type of genotype
# Checks if VCF, GWASpoly, or k-matrix.
#-------------------------------------------------------------
checkGenotypeFormat <- function (genotypeFile, ploidy) {
	message ("Checking genotype file format...")
	con = file (genotypeFile, "r")
	firstLine = readLines (con, n=1)
	close (con)

	# Check if VCF
	if (grepl ("VCF", firstLine)) {
		genotypeFile = convertVCFToACGTByNGSEP (genotypeFile) #output: filename.csv
		return (genotypeFile)
	}

	# Compare if all cells have the same lengths (nchars)
	# if True --> matrix type, False --> GWASpoly type
	equalLength <- function (cell, ploidy) {
		cellChar = as.character (cell)
		return (nchar (cellChar)==ploidy)
	}	

	# Check if k-matrix
	data   = read.csv (genotypeFile, row.names=1, check.names=F)
	sample = data [1:10,1:10]

	if (all(sapply(sample, equalLength, ploidy), na.rm=T)==TRUE) {
		N = nrow (data)
		newFilename = paste0 (strsplit (genotypeFile, "[.]")[[1]][1], "_GWASpoly.csv")
		genoGwaspoly = data.frame (Markers=rownames (data), Chrom=1, Position=1:N, data)
		write.csv (genoGwaspoly, newFilename, row.names=F)
		return (newFilename)
	}

	# Check if GWASpoly
	sample = data [1:10,3:10]
	if (all(sapply(sample, equalLength, ploidy), na.rm=T)==TRUE)
		return (genotypeFile)

	return ("Unknown-genotype")
}

#-------------------------------------------------------------
#-------------------------------------------------------------
#convertGenotypeMatrixToGwaspoly <- function (genotypeFile) {
source ("gwas-preprocessing.R")


args = commandArgs (trailingOnly=T)
genotypeFile = args [1]
ploidy = args [2]
	
print (checkGenotypeFormat (genotypeFile, ploidy))
