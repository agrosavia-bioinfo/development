#!/usr/bin/Rscript

# Construct a table with SNP|CHR|POS|REF|ALT
# Uses two files: sccafold (for Ref Alt) and genome (for chr, pos)

# Remove duplicated SNPs
source ("lglib.R")

# get-MAP-info-genome-scaffolds.R potato_infinium_8303_map_context_DM_v3_superscaffolds.txt potato_8303SNPs_potato_dm_v4.03.gff3
args = commandArgs (trailingOnly=T)

#----------------------------------------------------------
#----------------------------------------------------------
main <- function () {
	scaffoldFile = args [1]
	genomeFile   = args [2]

	refAltFile = getRefAltInfo (scaffoldFile)
	mapFile    = getSnpChromosomePositionFromGenome (genomeFile)

	refAlt  = read.csv (refAltFile)
	print (dim(refAlt))
	hd (refAlt)
	map     = read.csv (mapFile)
	print (dim(map))
	hd (map[])

	mapRefAlt = merge (x=map, y=refAlt, by="SNP", all.x=T)
	outFile   = paste0 (strsplit (scaffoldFile, "[.]")[[1]][1], "-MAP-REFALT.csv")
	write.csv (mapRefAlt, file=outFile, quote=F, row.names=F)
}


#----------------------------------------------------------
#----------------------------------------------------------
getRefAltInfo <- function (inFile) { 
	snpsAll   = read.table (inFile, header=T)
	snps      = snpsAll [!duplicated (snpsAll$Name),]
	ids       = snps [1]
	seq       = snps [4]
	alleles   = unlist(strsplit (substr (seq[,1], 52,54), split="/"))
	allelesDF = data.frame (matrix (alleles, ncol=2, byrow=T))
	ref      = as.character (allelesDF [,1])
	alt      = as.character (allelesDF [,2])

	out = cbind (SNP=ids, REF=ref, ALT=alt)
	colnames (out) = c("SNP", "REF", "ALT")
	refAltFile = paste0 (strsplit (inFile, split="\\.")[[1]][1], "-REFALT.csv")
	write.csv (file=refAltFile, out, quote=F, row.names=F)
	return (refAltFile)
}


#---------------------------------------------------------
# Get from Genome information file (.gff3) the SNP, CHROM, POS
# Write results to table "genome-SNPs-Chrom-Pos.tbl"
#---------------------------------------------------------
getSnpChromosomePositionFromGenome <- function (genomeFile) {
	#>>>>>>>>>>> local function <<<<<<<<<<<<<<<
	# Get SNP id from text string
	getIdMarker <- function (s) {
		return (strsplit (strsplit (s, split="=")[[1]][3], split=";")[[1]][1])
	}
	# Remove prefix "chr" from "chrXX"
	getChromMarker <- function (s) {
		return (strsplit (s, split="chr")[[1]][2])
	}
	#>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<

	# Get important columns from genome: CHRM, POS, SNP Info
	print (genomeFile)
	print (genomeFile)
	genomeAll    = read.table (file=genomeFile, sep="\t", header=F)
	genomeTmp = genomeAll [,c(1, 4, 9)]

	# Parsing of markers ids and chromosomes value
	genomeMarkers     = unlist (lapply (as.vector (genomeAll [,c(9)]), getIdMarker))
	#genomeMarkers     = gsub ("solcap_snp_","", genomeMarkers)  # To remove long prefix
	genomeChromosomes = as.numeric (unlist (lapply (as.vector (genomeAll [,c(1)]), getChromMarker)))
	genomePositions   = unlist (genomeAll [,c(4)])

	# Create new genome and write
	genome       = data.frame (SNP=genomeMarkers, CHROM=genomeChromosomes, POS=genomePositions)

	#gf=genotypeFiltered = gc [!duplicated (gc[,1],),]
	genomeNoDups = genome [!duplicated (genome$SNP),]
	genomeNoDups = genomeNoDups [order (genomeNoDups$SNP),]
	mapFile = paste0 (strsplit (genomeFile, "[.]")[[1]][1], "-MAP.csv")
	write.csv (file=mapFile, genomeNoDups, quote=F, row.names=F)
	return (mapFile)
} 

main ()
