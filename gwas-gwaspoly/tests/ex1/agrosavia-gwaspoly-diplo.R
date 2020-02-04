#!/usr/bin/Rscript
library (GWASpoly)

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
data2 <- set.K(data)
params = NULL

# Used to include population structure covariates
#params <- set.params(n.PC=10)
#params <- set.params(fixed=c("K1","K2","K3","K4"),
#                     fixed.type=rep("numeric",4))
# GWAS execution
message (">>> Running GWASpoly...")
#data3 = GWASpoly(data2, models=c("general","additive","1-dom", "2-dom"),traits=c("Gota"),n.core=4)
data3 = GWASpoly(data2, models=c("general","additive"),traits=c("Gota"),n.core=4)
# QQ-plot Output
message ("Ploting results...")
pdf (file="plots-qq-gwaspoly.pdf")
  par(mfrow=c(2,3)) #specifies a 2 x 3 panel
  #models <- c("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
  models <- c("additive","general")
  for (i in 1:2) {
    qq.plot(data3,trait="Gota",model=models[i])
  }
dev.off()

# QTL Detection
#data4 = set.threshold (data3, method="Bonferroni",level=0.05, n.core=4)
data4 = set.threshold (data3, method="FDR",level=0.05, n.core=4)
get.QTL (data4)

# Manhattan plot Output
pdf (file="plots-manhattan-gwaspoly.pdf")
  par(mfrow=c(2,3)) #specifies a 1 x 3 panel
  for (i in 1:2) {
    manhattan.plot (data4, trait="Gota", model=models [i])
  }  
dev.off ()


