
source ("lglib01.R")

genoFile = "tpg2plantgenome2015080073-sup-0002.csv"
geno = read.csv (genoFile, check.names=F)
hd (geno)

cnames = colnames (geno)
validNames = cnames [!grepl ("-", cnames)]

newGeno = geno [,validNames]
print (dim(geno))

newFilename = paste0(strsplit(genoFile,"[.]")[[1]][1], "-validNames.csv")
write.csv (newGeno, newFilename, row.names=F, quote=F)
