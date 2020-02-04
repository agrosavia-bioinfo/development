#!/usr/bin/Rscript

# Convert gwaspoly genotype format to plink transposed format (tped, tfam)

args = commandArgs(trailingOnly = TRUE)

args = c("genotype-diplo.tbl","phenotype.tbl")
gwaspolyGenotypeFile  = args [1]
gwaspolyPhenotypeFile = args [2]

gwaspolyGenotype = read.csv (file=gwaspolyGenotypeFile, header=T)
gw = gwaspolyGenotype

iid     = gwaspolyGenotype [,1]
fid     = iid
chrm    = gwaspolyGenotype [,2]+1
dist    = 0
pos     = gwaspolyGenotype [3]
alleles = gwaspolyGenotype [,-c(1,2,3)]

alleles [alleles==0] = "A A"
alleles [alleles==1] = "A B"
alleles [alleles==2] = "B B"
alleles [is.na(alleles)] = "0 0"

plinkGenotypeTPED = cbind (chrm,iid,dist,pos,alleles)
pl = plinkGenotypeTPED

plinkTPEDFilename = paste0 (strsplit (gwaspolyGenotypeFile, split="[.]")[[1]][1], "-plink.tped")
write.table (file=plinkTPEDFilename, plinkGenotypeTPED, col.names=F, row.names=F, quote=F)

snpsNames = colnames (gwaspolyGenotype)[-c(1,2,3)]
iid =fid = snpsNames
plinkFamilyTFAM = cbind (iid,fid,0,0,0,-9)
pf = plinkFamilyTFAM

plinkTFAMFilename = paste0 (strsplit (gwaspolyGenotypeFile, split="[.]")[[1]][1], "-plink.tfam")
write.table (file=plinkTFAMFilename, plinkFamilyTFAM, col.names=F, row.names=F, quote=F)

## Format phenotype
gwaspolyPhenotype = read.csv (file=gwaspolyPhenotypeFile, header=T)
idNames = as.character (gwaspolyPhenotype [,1])
IID = sapply (strsplit (idNames, split="And_",fixed=T), function (x) x[[2]][1])
FID = IID
BLIGHT = gwaspolyPhenotype [,2]
plinkPheno = cbind (FID, IID, BLIGHT)

plinkPhenoFilename = paste0 (strsplit (gwaspolyPhenotypeFile, split="[.]")[[1]][1], "-plink.tbl")
write.table (file=plinkPhenoFilename, plinkPheno, col.names=T, row.names=F, quote=F)
