#!/usr/bin/Rscript

#----------------------------------------------------------
# Convert numeric (gwaspoly) genotype to ACGT using solcap ref/alt alleles
# Basic genotype: [Markers+Alleles]
#----------------------------------------------------------
numericToACGTFormatGenotype <- function (genotypeFile, SNPsFile) 
{
	genotype     <- read.csv (genotypeFile, header=T, check.names=F)
	SNPs         <- read.csv (SNPsFile, header=T, check.names=F)
	rownames (SNPs) <- SNPs [,1]
	alleles      <- genotype [,-c(2,3)]
	allelesACGT  <- numericToACGTFormatAlleles (alleles, SNPs)
	genoACGT     = data.frame (genotype [,1:3], allelesACGT [,-1])

	outFile = paste0 (strsplit (genotypeFile, split="[.]")[[1]][1],"-ACGT.csv")
	#msgmsg ("Writing ACGT genotype to ", outFile, "...")
	write.csv (genoACGT, outFile, row.names=F)
}

numericToACGTFormatAlleles <- function (alleles, SNPs) 
{
	setA <- function (allelesVec, refs, alts) {
		id  = allelesVec [1]
		gnt <- as.numeric (allelesVec [-1])
		ref = alts [id,2]
		alt = refs [id,2]
		gnt [gnt==4] = strrep(ref,4)
		gnt [gnt==3] = paste0 (strrep(alt,1),strrep(ref,3))
		gnt [gnt==2] = paste0 (strrep(alt,2),strrep(ref,2))
		gnt [gnt==1] = paste0 (strrep(alt,3),strrep(ref,1))
		gnt [gnt==0] = strrep(alt,4)
		return (gnt)
	}
	refs <- data.frame (SNPs [,c(1,3)])
	rownames (refs) <- SNPs [,1]
	alts <- data.frame (SNPs [,c(1,2)])
	rownames (alts) <- SNPs [,1]
	alls <- alleles

	allelesNUM <- t(apply (alleles,  1, setA, refs, alts ))
	colnames (allelesNUM) = colnames (alleles [-1])
	rownames (allelesNUM) = rownames (alleles)

	newAlleles <- data.frame (Markers=alleles [,1], allelesNUM)
	return (newAlleles)
}

#----------------------------------------------------------
# Main
#----------------------------------------------------------
#args = c ("agrosavia-genotype-checked-tetra-NUM.tbl", "potato_infinium_8303_map_context_DM_v3_superscaffolds-ref-alt-snps.tbl")
args = commandArgs(trailingOnly = TRUE)
#args = c ("agrosavia-genotype-NUM-CLEANED-MAP.csv", "agrosavia-genotype-MAP-altref.csv")

genotypeFile  = args [1]
mapSNPs  = args [2]

numericToACGTFormatGenotype (genotypeFile, mapSNPs)

