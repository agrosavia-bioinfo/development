#!/usr/bin/Rscript
#-----------------------------------------------------------
# Create genotype by checking Marker SNPs in the phenotype file
#-----------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

# inputs
inGenoFile = args [1]
inPhenoFile  = "agrosavia-phenotype-gwaspoly.tbl"

#inPhenoFile = "agrosavia-phenotype-gwaspoly.tbl"
#inGenoFile  = "agrosavia-genotype-gwaspoly.tbl"

# outputs
filename = strsplit (inGenoFile, "[.]")[[1]][1] 
outGenotypeGwaspoly = sprintf ("%s-checked.tbl", filename)
outInvalidSNPs      = sprintf ("%s-checked-errors.tbl", filename)


ph = read.table (inPhenoFile, header=T, sep=",")
snps = as.vector (ph [,1])

gn = read.table (inGenoFile, header=T, sep=",")
gnSnps = colnames (gn)[-(1:3)] 

gnTbl = data.frame (gn [,1:3])
errorsTbl = c ("Invalid SNPs: in Phenotype but not in Genotype")
for (sn in snps) {
  if (sn %in% gnSnps) {
    gnCol = data.frame (gn [, sn])
    colnames (gnCol) = sn
    gnTbl = cbind (gnTbl, gnCol)
  }else {
    print (sprintf (">>> Invalid Snp %s ", sn))
    errorsTbl = append (errorsTbl, sn)
  }
  
}

write.table (file=outGenotypeGwaspoly, gnTbl, quote=F, sep=",", row.names=T)
write.table (file=outInvalidSNPs, errorsTbl, col.names=F,row.names=T)
