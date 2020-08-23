#!/usr/bin/Rscript

library (vcfR)
source ("lglib01.R")

filename = "vt.vcf"
vcf = read.vcfR (filename)
head (vcf)

# Extract genotipe strings" 
gt = extract.gt (vcf)
head (gt)


# Extract map info (id, chrom, pos)
map = vcf@fix [,c(1:2,4:5)]
head (map)

# Merge gt + info
gtmap = cbind (SAMPLES=rownames (gt), map, gt)

hd(gtmap)

# ">>>>>>>>>>>"
vcfToACGT <- function (allelesMat, refsVec, altsVec) {
	allele = allelesMat [1]
	ref    = refsVec [1]
	alt    = altsVec [1]

	nums   = strsplit (allele, split=c("|", "\\\\")

	gnt [gnt==4] = strrep(ref,4)
	gnt [gnt==3] = paste0 (strrep(alt,1),strrep(ref,3))
	gnt [gnt==2] = paste0 (strrep(alt,2),strrep(ref,2))
	gnt [gnt==1] = paste0 (strrep(alt,3),strrep(ref,1))
	gnt [gnt==0] = strrep(alt,4)
	return (gnt)
}
# ">>>>>>>>>>>"
allelesMat = gtmap [,-1:-5]
refsVec    = gtmap [,4]
altsVec    = gtmap [,5]

vcfToACGT (allelesMat, refsVec, altsVec)






