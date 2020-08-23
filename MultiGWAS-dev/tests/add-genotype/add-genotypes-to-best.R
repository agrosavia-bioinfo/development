#!/usr/bin/Rscript

# Add genotypes to selected SNPs 

source ("lgRlib.R")

selFilename       = args [1]
genoFilenameACGT  = args [2]
genoFilenameNUM   = args [3]


sel = read.csv (file=selFilename, sep="\t", stringsAsFactors=F)
hd (sel)

genoACGT = read.csv (file=genoFilenameACGT)
hd (genoACGT)
print (dim(genoACGT))

genoNUM = read.csv (file=genoFilenameNUM)
hd (genoNUM)
print (dim(genoNUM))


# for ACGT
geno = genoACGT
alleles = geno [match (sel$SNP, geno$Marker),]; 
hd (alleles,,n=13)

selgeno = cbind (sel, alleles[,-c(1:3)])
hd (selgeno)

outFilename = paste0 (strsplit (selFilename, "[.]")[[1]][1], "-withGenotypes-ACGT.csv")
write.table (file=outFilename, selgeno, sep="\t", quote=F, row.names=F)

# for NUM
geno = genoNUM
alleles = geno [match (sel$SNP, geno$Marker),]; 
hd (alleles,,n=13)

selgeno = cbind (sel, alleles[,-c(1:3)])
hd (selgeno)

outFilename = paste0 (strsplit (selFilename, "[.]")[[1]][1], "-withGenotypes-NUM.csv")
write.table (file=outFilename, selgeno, sep="\t", quote=F, row.names=F)


