#!/usr/bin/Rscript

# Extract  the common samples from genotype and phenotypes
# and write to new files

source ("lglib.R")

#-------------------------------------------------------------
# Add label to filename and new extension (optional)
#-------------------------------------------------------------
addLabel <- function (filename, label, newExt=NULL)  {
	nameext = strsplit (filename, split="[.]")
	name    = nameext [[1]][1] 
	if (is.null (newExt))
		ext     = nameext [[1]][2] 
	else
		ext     = newExt
	newName = paste0 (nameext [[1]][1], "-", label, ".", ext )
	return (newName)
}
filterByCommonNames <- function (genotypeFile, phenotypeFile) 
{
	geno  = read.csv (file=genotypeFile, header=T)
	pheno = read.csv (file=phenotypeFile, header=T)
	rownames (pheno) = pheno [,1]

	genoColumns   = colnames (geno)
	phenoSamples  = pheno [,1] 
	commonSamples <- intersect (genoColumns, phenoSamples) 

	genoCommon  <- geno  [,c(genoColumns[1:3], commonSamples)]
	hd (genoCommon)
	map = genoCommon [, (1:3)]
	write.table (file="map.tbl", map, quote=F, row.names=F, sep="\t")

	print (length (commonSamples))


	phenoCommon <- pheno [pheno[,1] %in% commonSamples,]
	hd (phenoCommon)
	trait  <- colnames (phenoCommon)[2]
	genoCommonFile  = addLabel (genotypeFile, "COMMON")
	phenoCommonFile = addLabel (phenotypeFile, "COMMON")
	write.table (file=genoCommonFile, genoCommon, quote=F, row.names=F, sep=",")
	write.table (file=phenoCommonFile, phenoCommon, quote=F, row.names=F, sep=",")

	return (list (genotypeFile=genoCommonFile, phenotypeFile=phenoCommonFile, trait=trait))
}



filterByCommonNames (args [1], args [2])
