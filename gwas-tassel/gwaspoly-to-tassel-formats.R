#!/usr/bin/Rscript

# Convert gwaspoly genotype format to plink transposed format (tped, tfam)

#------------------------------------------------------------------------------
## Write and format tassel structure
#------------------------------------------------------------------------------
writeTasselStruct <- function (structureFile) {
	structData = read.csv (file=structureFile, header=T)
	idNames = as.character (structData [,1])
	Trait = sapply (strsplit (idNames, split="And_",fixed=T), function (x) x[[2]][1])
	struct = structData [,-1]
	newData = cbind ("<Trait>"=Trait, struct)

	message (">>> Writing tassel structure")
	outFilename = paste0 (strsplit (structureFile, split="[.]")[[1]][1], "-tassel.tbl")
	sink (outFilename)
	cat ("<Covariate>\n")
	write.table (file="", newData, col.names=T, row.names=F, quote=F, sep="\t")
	sink()
}

#------------------------------------------------------------------------------
## Format and write tassel phenotype
#------------------------------------------------------------------------------
writeTasselPheno <- function (gwaspolyPhenotypeFile) {
	gwaspolyPhenotype = read.csv (file=gwaspolyPhenotypeFile, header=T)
	idNames = as.character (gwaspolyPhenotype [,1])
	Taxa = sapply (strsplit (idNames, split="And_",fixed=T), function (x) x[[2]][1])
	BLIGHT = gwaspolyPhenotype [,2]
	plinkPheno = cbind (Taxa, BLIGHT)

	message (">>> Writing tassel phenotype...")
	plinkPhenoFilename = paste0 (strsplit (gwaspolyPhenotypeFile, split="[.]")[[1]][1], "-tassel.tbl")
	sink (plinkPhenoFilename)
	cat ("<Phenotype>\n")
	cat ("taxa\tdata\n")
	write.table (file="", plinkPheno, col.names=T, row.names=F, quote=F, sep="\t")
	sink()
}

#----------------------------------------------------------
# Transform a numeric genotype (0,1,2) to ACGT format
#----------------------------------------------------------
transformGenotypeABtoGC <- function (genotypeFilename, snpsRefAltFile) {
	changeABtoGC <- function (snp, ref, alt) {
		alleles = snp [-c(1,2,3)]
		alleles [alleles==0] = paste0(ref,ref)
		alleles [alleles==1] = paste0(ref,alt)
		alleles [alleles==2] = paste0(alt,alt)
		alleles [is.na(alleles)] = "00"
		snp [-c(1,2,3)] = alleles
		return (snp)

	}

	snpsRefAlt <- read.table (file=snpsRefAltFile, header=T)
	snpsRefAlt <- snpsRefAlt [!duplicated (snpsRefAlt [,1]),]
	genotype <- read.csv (file=genotypeFilename)

	print (nrow (genotype))
	for (i in 1:(nrow (genotype))) {
		id  = as.character (genotype [i,1])
		ref = snpsRefAlt [snpsRefAlt==id,2]  
		alt = snpsRefAlt [snpsRefAlt==id,3]  
		genotype [i,] = changeABtoGC (genotype [i,], ref, alt)
		#genotype [i,] = gi
	}
	message (">>> Writing new GC genotype...")
	name = strsplit (genotypeFilename, split="\\.")[[1]][1]
	newGenotypeFilename = paste0 (name, "-AG.tbl")
	write.table (file=newGenotypeFilename, genotype,  quote=F)
	return (newGenotypeFilename)
}

#----------------------------------------------------------
# Add tabs to alleels changign AG --> A	G
#----------------------------------------------------------
transformAlleles <- function (alleles) {
	ncols = ncol (alleles)
	nrows = nrow (alleles)
	alleles = t (alleles)
	tabs <- function (x) {return (sprintf("%s\t%s", substr(x,1,1),substr(x,2,2)))}
	allelesTabbed =matrix (sapply (alleles,tabs), ncol=nrows, nrow=ncols, byrow=T)
	return (allelesTabbed)
}

#----------------------------------------------------------
# Main
#----------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)

args = c("genotype-diplo.tbl","phenotype.tbl", "potato_infinium_8303_map_context_DM_v3_superscaffolds-ref-alt-snps.tbl",
		 "structure-checked.tbl")
gwaspolyGenotypeFile  = args [1]
gwaspolyPhenotypeFile = args [2]
snpsRefAltFile        = args [3]
structureFile         = args [4]

#genotypeFilename = transformGenotypeABtoGC (gwaspolyGenotypeFile, snpsRefAltFile)
genotypeFilename = "genotype-diplo-AG.tbl"

gwaspolyGenotype = read.table (file=genotypeFilename, header=T)
gw = gwaspolyGenotype

snpsIds = gsub ("solcap_snp_", "", gwaspolyGenotype [,1])
fid     = snpsIds
chrm    = gwaspolyGenotype [,2]+1
dist    = 0
pos     = gwaspolyGenotype [3]
alleles = gwaspolyGenotype [,-c(1,2,3)]
samplesIds = gsub ("And_", "", colnames (alleles))

message (">>> Writing PED files...")
talleles = transformAlleles (alleles)
namesGeno = c("fid", "iid", "pid", "mid", "sex", "phe", rep ("X", ncol (talleles)))
plinkGenotypePED = cbind (-9, samplesIds, -9, -9, -9, -9, talleles)
colnames (plinkGenotypePED) = namesGeno
plinkPEDFilename = paste0 (strsplit (gwaspolyGenotypeFile, split="[.]")[[1]][1], "-plink.ped")
write.table (file=plinkPEDFilename, plinkGenotypePED, col.names=F, row.names=F, quote=F, sep="\t")

message (">>> Writing MAP files...")
plinkFamilyMAP = cbind (chr=chrm, iid=snpsIds, dist=-9, pos=pos)
plinkMAPFilename = paste0 (strsplit (gwaspolyGenotypeFile, split="[.]")[[1]][1], "-plink.map")
write.table (file=plinkMAPFilename, plinkFamilyMAP, col.names=F, row.names=F, quote=F, sep="\t")

## Write and format tassel phenotype
writeTasselPheno (gwaspolyPhenotypeFile)

## Write and format tassel structure
writeTasselStruct (structureFile) {


