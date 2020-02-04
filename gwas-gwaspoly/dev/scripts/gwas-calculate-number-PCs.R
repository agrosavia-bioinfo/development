#!/usr/bin/Rscript

# Calculate the number or PCs from Potato Dossage Matrix
# - Load the genotype data and 
# - Extract the dosage data 
# - Calculate the max number of PCs
  
library (missMDA) # For number of PCA with NAs
library (FactoMineR) # For PCA

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) {
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}

runPCA <- function (data) {
  outPCA = PCA (data)
  
  print (outPCA$sdev[1:10])
  plot (outPCA)
  
  return (outPCA)
}

runPCAWithNAs <- function (data) {
  #Estimate the number of dimensions for the PCA by cross-validation
  nb = estim_ncpPCA (data, ncp.max=5)
  res.comp = imputePCA(data, ncp=2)
  res.pca  = PCA (res.comp)
  
  print (res.comp$sdev[1:10])
  plot (res.comp, type="l")
}

inGenotypeFilename     = "agrosavia-genotype-gwaspoly-checked.tbl"
genotype               = read.table (inGenotypeFilename, header=T, sep=",")
dossage                = genotype [,-(1:3)]
dossage.noNAs          = na.omit (dossage)
dossage.transposed     = t (dossage)
tpd                    = dossage.transposed

tpd.noNAs              = na.omit (tpd)
print (dim(tpd.noNAs))


dsPCA = runPCA (dossage.noNAs)

msg ("Dimensions of dossage:", dim(dossage))
msg ("Dimensions of dossage transposed:", dim (dossage.transposed))

dst = runPCAWithNAs(dossage)
