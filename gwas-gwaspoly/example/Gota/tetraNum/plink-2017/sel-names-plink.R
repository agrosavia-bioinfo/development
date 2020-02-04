#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE)
args = c("agrosavia-genotype-FLT-plink-bin.bim", "agrosavia-genotype-FLT-plink-bin.fam")
markersFile = args [1]
samplesFile = args [2]
 
md = read.table (file=markersFile)
sd = read.table (file=samplesFile )

markers = md [,2]
samples = paste0 ("A", sd [,2])

write.table (file="markers2017.tbl", markers, row.names=F, quote=F, col.names=F)
write.table (file="samples2017.tbl", samples, row.names=F, quote=F, col.names=F)

# Select from geno/pheno
genoFile  = "agrosavia-genotype-checked-tetra-NUM.tbl"
phenoFile = "agrosavia-phenotype-checked-Gota.tbl"

geno  = read.csv (file=genoFile, header=T)
rownames (geno) = geno [,1]
pheno = read.csv (phenoFile, header=T)
rownames (pheno) = pheno [,1]
head (pheno)
head (samples)

newGeno = geno [markers,]
newGenoName = "agrosavia-genotype-2017.tbl"
write.csv (file=newGenoName, newGeno, quote=F, row.names=F)

newPheno = pheno [samples,]
head (newPheno)
newPhenoName = "agrosavia-phenotype-2017.tbl"
write.csv (file=newPhenoName, newPheno, quote=F, row.names=F)



