#!/usr/bin/Rscript
library (GWASpoly)

TRAIT="tuber_shape"

snpsModelsDiplo = c("general","additive","1-dom")
testModelsDiplo = c("general", "additive","1-dom-alt","1-dom-ref")
nPlots          = 4
#-------------------------------------------------------------
# Add label to filename
#-------------------------------------------------------------
addLabel <- function (filename, label)  {
	nameext = strsplit (filename, split="[.]")
	newName = paste0 (nameext [[1]][1], "-", label, ".", nameext [[1]][2])
	return (newName)
}

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) {
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}
#-------------------------------------------------------------
  
args = commandArgs(trailingOnly = TRUE)
print (args)

genotypeFile  = args [1]
phenotypeFile = args [2]

ph = read.table (phenotypeFile, header=T, sep=",")
gn = read.table (genotypeFile, header=T, sep=",")


# Read input genotype and genotype (format: "numeric" or "ACGT")
message (">>> Reading data...")
data = read.GWASpoly (ploidy = 2, delim=",", format = "ACGT", n.traits = 1, 
                      pheno.file = phenotypeFile, geno.file = genotypeFile)

message (">>> Calculating kinship...")
# Populations structure by kinship
data2 <- set.K(data, K=NULL)

# Used to include population structure covariates
#params <- set.params(n.PC=10)
#params <- set.params(fixed=c("K1","K2","K3","K4"),
#                     fixed.type=rep("numeric",4))
params = NULL
# GWAS execution
message (">>> Running GWASpoly...")
data3 = GWASpoly(data2, models=snpsModelsDiplo,traits=c(TRAIT),n.core=4)


# QTL Detection
#data4 = set.threshold (data3, method="Bonferroni",level=0.05, n.core=4)
data4 = set.threshold (data3, method="FDR",level=0.05, n.core=4)
significativeQTLs = get.QTL (data4)
outFile = paste0 ("out-QTL-", TRAIT,".scores")
write.table (file=outFile, significativeQTLs, quote=F, sep="\t", row.names=F)

# QQ-plot Output
message ("Ploting results...")
outFile = paste0 ("out-QQ-", TRAIT,".pdf")
pdf (file=outFile)
  par(mfrow=c(2,3)) #specifies a 2 x 3 panel
  models <- c(testModelsDiplo)
  for (i in 1:nPlots) {
    qq.plot(data3,trait=TRAIT,model=models[i])
  }
dev.off()


# Manhattan plot Output
outFile = paste0 ("out-manhattan-", TRAIT,".pdf")
pdf (file=outFile)
  par(mfrow=c(2,3)) #specifies a 1 x 3 panel
  for (i in 1:nPlots) {
    manhattan.plot (data4, trait=TRAIT, model=models [i])
  }  
dev.off ()



