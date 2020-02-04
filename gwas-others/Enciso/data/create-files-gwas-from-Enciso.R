#!/usr/bin/Rscript

# Create phenotypes by year processing the original phenotype


args = commandArgs(trailingOnly = TRUE)
options (width=300)

args = c ("genotype-enciso.tbl", "phenotype-enciso.tbl", "genome.gff3")
genotypeFile   = args [1]
phenotypeFile  = args [2]
genomeFile     = args [3]

#-------------------------------------------------------------
# Shorcut util functions
#-------------------------------------------------------------
s <- function (s,cols=3)
	print (s[1:5,1:cols])

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) {
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}

#----------------------------------------------------------
# Create phenotypes for year 
#----------------------------------------------------------
createPhenotypeByYear <- function (phenotypeFile){
	phenotypeAll  = read.table (file=phenotypeFile, header=T)
	for (year in c(2010:2015, 2017)) {
		print (year)
		phenoByYearAll = phenotypeAll [phenotypeAll$Year==year,]
		phenoByYear    = aggregate (phenoByYearAll$Score, list (phenoByYearAll$Genotype), mean)
		phenoByYear    = aggregate (phenoByYearAll$Score, list (phenoByYearAll$Genotype), mean)
		colnames (phenoByYear) = c ("Stub", "LBlight")

		outfile = sprintf ("phenotype-Enciso-LBlight-%s.tbl", year)
		write.table (phenoByYear, file=outfile, quote=F, row.names=F, sep=",")
	}
}

#----------------------------------------------------------
# Extract genome annotations 
#----------------------------------------------------------
extractGenomeAnnotations <- function (genomeFile) {
	genome    = as.matrix (read.table (file=genomeFile, header=F))
	n         = nrow (genome)
	genomeAnn = matrix (ncol=3, nrow=n) 
	for (i in 1:nrow (genome)) {
		snp       = genome [i,]
		chr       = gsub ("chr(.+)", "\\1", snp[1])
		pos       = snp [4]
		name      = gsub (".*Name=(.+);.*", "\\1", snp [9])
		genomeAnn [i,] = c (name, chr, pos)
	}
	colnames (genomeAnn) = c("Marker", "Chrom", "Position")
	#write.table (file="genome-markers.tbl", genomeAnn, quote=F, sep="\t", row.names=F)
	return (data.frame (genomeAnn, row.names=NULL))
}

#----------------------------------------------------------
# Create GWASpoly genotype
#----------------------------------------------------------
createGWASpolyGenotype <- function  (genotypeFile, genomeAnn) {
	genotypeAll    = read.table (file=genotypeFile, header=T)
	genotypeBySNPs = t(genotypeAll)
	colnames (genotypeBySNPs) = rownames (genotypeAll)
	rownames (genotypeBySNPs) = colnames (genotypeAll)
	write.csv (file="genotype-snps-by-samples.tbl", genotypeBySNPs, row.names=T)

	genotypeBySNPsDframe = data.frame (Marker=rownames (genotypeBySNPs), genotypeBySNPs, row.names=NULL, check.names=F)
	gdf = genotypeBySNPsDframe
	genotypeChromPos = merge (genomeAnn, genotypeBySNPsDframe, by.x="Marker", by.y="Marker")
	gc = genotypeChromPos
	write.csv (file="genotype-Enciso-gwaspoly.tbl", genotypeChromPos, quote=F, row.names=F)
	return (genotypeChromPos)
}

#----------------------------------------------------------
#----------------------------------------------------------

msg ("Extracting annotations from genome...")
genomeAnn = extractGenomeAnnotations (genomeFile)

msg ("Creating GWASpoly genotype...")
genotype  = createGWASpolyGenotype (genotypeFile, genomeAnn)

createPhenotypeByYear (phenotypeFile)
