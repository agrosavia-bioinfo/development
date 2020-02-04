#!/usr/bin/Rscript

# Create the barplots for SNPs from high and low blight resistance
# It takes info from phenotype, genotype, and QTLs

library (dplyr)
options (stringAsFactors=F, width=300)

#----------------------------------------------------------
# Giving the type of phenotypes (low or high) it calculates the
# allele frequencies and plot them
#----------------------------------------------------------
conteo <- function (genotype, phenotypes, labelPheno, model, snp,i) {

	alleles = as.numeric ( (genotype %>% filter (Markers==snp) %>% as.data.frame())[,phenotypes])
	counts = (sapply (c(0,1,2,3,4), function (x,n) {length (which (x==n))}, alleles))

	alleleNames = c("AAAA","AAAB","AABB", "ABBB", "BBBB")
	df = as.data.frame (counts, row.names=alleleNames)
	name = sprintf("MODEL-%s-SNP-%s", model,snp)
	#pdf (sprintf("%sResistance-%s.pdf", labelPheno, name))
	barplot(df[,1],names=alleleNames, main=sprintf ("%s-%s",labelPheno, name), col=rep(i,5))
	#dev.off()
}

#----------------------------------------------------------
# Read the genotype and select two groups of 50 samples (high and low)
#----------------------------------------------------------
# Read and plot phenotype
selectHighLowPhenotypes <- function (phenotypeFile) {
	message (">>> Reading Genotype...")
	pheno = read.csv ("agrosavia-phenotype.tbl", header=T)
	ph <- pheno

	highSamples = data.frame (ph %>% arrange (desc(gota)) %>% head (50))
	hs <- as.data.frame (highSamples [,2])
	lowSamples = data.frame (ph %>% arrange (desc(gota)) %>% tail (50))
	ls <- as.data.frame (lowSamples [,2])

	pdf ("plot-lblight-resistance.pdf")
	boxplot (as.data.frame(c(hs,ls),col.names=c("High", "Low")), main="Late Blight Resistance" )
	#boxplot (as.data.frame(c(hs,ls),col.names=c("High Resistance", "Low Resistance")))
	dev.off()

	return (list (high=as.character(highSamples[,1]), low=as.character(lowSamples[,1])))
}
#----------------------------------------------------------
# Get most significative SNPs (QTLs)
# Return dataframe with Models and Markers
#----------------------------------------------------------
getSignificativeQTLs <- function (qtlsFile) {
	qtlsAll = read.table (qtlsFile, header=T, sep="\t")
	qtls = as.data.frame (qtlsAll %>% filter (GC >= 1.0))
	bestQTLs <- qtls %>% select (Model, Marker) %>% as.data.frame()

	return (bestQTLs)
}

#----------------------------------------------------------
# Create a matrix of barplots for each significative QTL
#----------------------------------------------------------
plotQTLsAlleles <- function (outFile, genotype, phenotypes, bestQTLs) {
	pdf (outFile, width=22, height=7)
	op=par (mfrow = c(2, nrow(bestQTLs)))
	colors = c(1,2,3,4,5,6)
	for (i in 1:nrow (bestQTLs)) {
		model = as.character (bestQTLs[i,1])
		snp = as.character (bestQTLs[i,2])
		conteo (genotype, phenotypes$high, "High", model, snp,i)
	}
	for (i in 1:nrow (bestQTLs)) {
		model = as.character (bestQTLs[i,1])
		snp = as.character (bestQTLs[i,2])
		conteo (genotype, phenotypes$low, "Low", model, snp,i)
	}
	par(op)
	dev.off()
}
#----------------------------------------------------------
# Main
#----------------------------------------------------------
phenotypeFile = "agrosavia-phenotype.tbl"
genotypeFile  = "agrosavia-genotype-renamed.tbl"
qtlsFile      = "qtls.tbl"
outFile       = "plot-alleles-significative-QTLs.pdf"

message (">>> Reading Phenotype")
phenotypes   = selectHighLowPhenotypes (phenotypeFile)

message (">>> Reading Genotype")
genotype     = read.csv (genotypeFile, header=T)

message (">>> Reading QTLs")
bestQTLs = getSignificativeQTLs (qtlsFile)

message (">>> Plotting Alleles from QTLs")
plotQTLsAlleles (outFile, genotype, phenotypes, bestQTLs)

