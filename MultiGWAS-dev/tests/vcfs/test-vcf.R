#!/usr/bin/Rscript


convertVCFtoACGT <- function (filename, outFilename="") {
	suppressMessages (library (vcfR))
	vcf = read.vcfR (filename, verbose=F)

	# Extract genotipe strings" 
	gt = extract.gt (vcf, return.alleles=T)

	# Extract map info (id, chrom, pos) Merge gt + info
	map = vcf@fix [,c(1:2,4:5)]
	gtmap = cbind (MARKERS=rownames (gt), map, gt)

	# >>>>>>>> Convert VCF Genotypes to ACGT. First by row, then by cell
	vcfToACGTForCell <- function (allelesCell) {
		if (is.na (allelesCell))
			return (NA)

		numsall    = strsplit (allelesCell, split="[|/]")
		str = sapply (numsall[[1]], substring, 1, 1)
		#genotype = sapply (
		#ref  = substring (nums[[1]][1], 1, 1)
		#alt  = substring (nums[[1]][2], 1, 1)
		#return (paste0(ref,alt))
		return (paste (str, collapse=""))
	}
	vcfToACGTForRow <- function (allelesRow) {
		rows = sapply (allelesRow, vcfToACGTForCell)
		return (rows)
	} # ">>>>>>>>>>>"

	allelesMat = gtmap [,-1:-5]
	allelesACGT <- t(apply (allelesMat, 1, vcfToACGTForRow))
	colnames (allelesACGT) = colnames (allelesMat)
	rownames (allelesACGT) = rownames (allelesMat)

	newAlleles <- data.frame (gtmap[,1:3], allelesACGT)
	if (outFilename=="")
		outFilename = paste0 (strsplit (filename, "[.]")[[1]][1], ".csv")

	write.csv (newAlleles, outFilename, row.names=F, quote=F)
	return (outFilename)
}


args = commandArgs(trailingOnly = TRUE)
filename = args [1]
convertVCFtoACGT (filename, "t.csv")



