
# Sort alleles from genotype
source ("lglib01.R")
filename = "ClusterCall_prediction_CCC.csv"
geno = read.csv (filename, check.names=F); hd (geno)
alleles = geno [,-1]; hd (alleles)
sortedNames = sort (colnames(alleles));
sortedAlleles = alleles [, sortedNames];hd (sortedAlleles)


genoSorted = cbind (Markers=geno [,1], sortedAlleles); hd (genoSorted)
outfile = addLabel (filename, "SORTED")
write.csv (genoSorted, outfile, row.names=F, quote=F)
