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
#args = c (phenotypeFile = "agrosavia-phenotype-gwaspoly.tbl", genotypeFile  = "agrosavia-genotype-gwaspoly-checked.tbl")
args = c ("phenotype-checked.tbl", "genotype-checked.tbl")
print (args)
phenotypeFile = args [1]
genotypeFile  = args [2]

testModels <- c ("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
#testModels  = c ("additive")
testTraits <- c ("gota")
nTraits=1
#testTraits  = c ("tuber_shape")


# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
message (">>> Reading data...")
data = read.GWASpoly (ploidy = 4, pheno.file = phenotypeFile, geno.file = genotypeFile, format = "AB", n.traits = nTraits, delim=",")

# Populations structure by kinship
message (">>> Calculating kinship...")
data2 <- set.K(data)
#params <- set.params(fixed=c("Grp1","Grp2","Grp3","Grp4"), fixed.type=rep("numeric",4))
params = NULL
#params <- set.params(n.PC=10)

# GWAS execution
message (">>> Running GWASpoly...")
data3 = GWASpoly(data2, models=c("general","additive","1-dom", "2-dom"),traits=c("gota"), params=params)

message (">>> Plotting results...")
# QQ-plot Output
pdf (file="plots-qq-gwaspoly.pdf")
  par(mfrow=c(2,3)) #specifies a 2 x 3 panel
  models <- testModels
  for (i in 1:length(testModels)) {
    qq.plot(data3,trait=testTraits [1], model=models[i])
  }
dev.off()

# QTL Detection
data4 = set.threshold (data3, method="Bonferroni",level=0.05)
get.QTL (data4)

# Manhattan plot Output
pdf (file="plots-manhattan-gwaspoly.pdf")
  par(mfrow=c(2,3)) #specifies a 1 x 3 panel
  models <- testModels
  for (i in 1:length(testModels)) {
    manhattan.plot (data4, trait=testTraits[1], model=models [i])
  }  
dev.off ()


