#!/usr/bin/Rscript
## Get a summarize table from GWASpoly results
options (width=300)
library(dplyr)

# Read data
nm=naiveModel     = read.table (file="naive.tbl", header=T, sep="\t")
km=kinshipModel   = read.table (file="kinship.tbl", header=T, sep="\t")
pm=structModel    = read.table (file="pcs.tbl", header=T, sep="\t")
kpm=kinStructModel = read.table (file="kinship+pcs.tbl", header=T, sep="\t")
#qm=structModel    = read.table (file="structure.tbl", header=T, sep="\t")
#kqm=kinStructModel = read.table (file="kinship+structure.tbl", header=T, sep="\t")

# Organize data
nmt  = cbind (Ploidy="4", GWAS="Naive",  as.data.frame (nm  %>% group_by (Model) %>% mutate (absGC = abs(GC)-1) %>% arrange (absGC) %>% summarize (GC=mean(GC), nSNPs=n())))
kmt  = cbind (Ploidy="4", GWAS="Kinship",  as.data.frame (km  %>% group_by (Model) %>% mutate (absGC = abs(GC)-1) %>% arrange (absGC) %>% summarize (GC=mean(GC), nSNPs=n())))
pmt  = cbind (Ploidy="4", GWAS="PCs",  as.data.frame (pm  %>% group_by (Model) %>% mutate (absGC = abs(GC)-1) %>% arrange (absGC) %>% summarize (GC=mean(GC), nSNPs=n())))
kpmt = cbind (Ploidy="4", GWAS="Kin+PCs", as.data.frame (kpm %>% group_by (Model) %>% mutate (absGC = abs(GC)-1) %>% arrange (absGC) %>% summarize (GC=mean(GC), nSNPs=n())))
#qmt  = cbind (Ploidy="4", GWAS="Structure",  as.data.frame (qm  %>% group_by (Model) %>% mutate (absGC = abs(GC)-1) %>% arrange (absGC) %>% summarize (GC=mean(GC), nSNPs=n())))
#kqmt = cbind (Ploidy="4", GWAS="Kin+Str", as.data.frame (kqm %>% group_by (Model) %>% mutate (absGC = abs(GC)-1) %>% arrange (absGC) %>% summarize (GC=mean(GC), nSNPs=n())))

print (">>> nmt:")
print (nmt)

smmTable  = rbind (nmt, kmt, pmt, kpmt) 
#smmTable  = rbind (nmt, kmt, qmt, kqmt) 
message ("")
print (">>> smm:")
print (smmTable)

# Write data summarization
write.table (smmTable, "gwaspoly-summary-models.tbl", quote=F, sep="\t", row.names=F)
