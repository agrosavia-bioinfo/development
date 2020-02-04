#!/usr/bin/Rscript


suppressMessages (library (dplyr))
options (width=300, stringsAsFactors=T)
args = commandArgs(trailingOnly=T)


#----------------------------------------------------------
# Plink separated alleles (e.g. A C A G A A) are joined (AC AG AA)
#----------------------------------------------------------
spaceAlleles <- function (alleles) {
	all <<-alleles
	sep <- function (a) {
		s="";
		for (i in 1:4) 
			s=paste0 (s, substr (a,start=i,stop=i)," ");
		return (s)
	}

	allelesSep <<- sapply (as.character(alleles), sep)
	allelesMat <<- matrix (allelesSep, nrow=nrow(alleles), ncol=ncol(alleles), byrow=F)
	return (allelesMat)
}


args = c("agrosavia-genotype-tetra-ACGT.tbl", "agrosavia-phenotype.tbl", 
		 "samples-filtered-names.tbl", "markers-filtered-names.tbl")
genoFile    = args [1]
phenoFile   = args [2]
samplesFile = args [3]
markesFile  = args [4]


geno    = read.table (file=genoFile, header=T, row.names=1)
pheno   = read.table (file=phenoFile, header=T, row.names=1, sep=",")
samples = read.table (file=samplesFile, header=T)
samples = as.character (samples[,1])
markers = read.table (file=markesFile, header=T)
markers = as.character (markers[,1])



phenoFiltered = cbind (Samples=samples, Trait=pheno [samples,])
write.table (file="phenoFiltered.tbl", phenoFiltered, quote=F,row.names=F, sep="\t")

genoFilteredMarkers  = cbind (Markers=markers, geno[markers,c(samples)])# %>% filter(geno$Markers %in% markers)
write.table (file="genoFilteredByMarkers.tbl", genoFilteredMarkers, quote=F,row.names=F, sep="\t")

genoFilteredSamples = cbind (Sample=samples, Trait=phenoFiltered[,2],  t(geno[markers,c(samples)]))
write.table (file="genoFilteredBySamples.tbl", genoFilteredSamples, quote=F,row.names=F, sep="\t")

alleles = geno[markers,c(samples)]
spacedAlleles = spaceAlleles (t(alleles))
genoFilteredSamples = cbind (Sample=samples, Trait=phenoFiltered[,2],  spacedAlleles)
write.table (file="genoFiltered-SHEsis.tbl", genoFilteredSamples, quote=F,row.names=F,col.names=F, sep="\t")

#genoFilteredSamples = cbind (Sample=samples, Trait=phenoFiltered[,2], spacedAlleles)
#genoFiltered  = cbind (samples, t (genoFiltered))



