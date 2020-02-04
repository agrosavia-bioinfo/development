#!/usr/bin/Rscript

#----------------------------------------------------------
# Execute GWASPoly with example pheno and geno files
#----------------------------------------------------------

library (GWASpoly)
TRAIT="Gota"

args = commandArgs(trailingOnly = TRUE)
#args = c ("agrosavia-genotype-checked.tbl", "agrosavia-phenotype-checked.tbl")
args = c ("agrosavia-genotype-tetra-ACGT.tbl", "agrosavia-phenotype-checked.tbl")
print (args)
genotypeFile  = args [1]
phenotypeFile = args [2]

data = read.GWASpoly (ploidy = 4, pheno.file = phenotypeFile, geno.file = genotypeFile, format = "ACGT", n.traits = 1, delim=",")

#ph = read.table (phenotypeFile, header=T, row.names = 1, sep=",")
#gn = read.table (genotypeFile, header=T, row.names = 1, sep=",")

##------------GWAS Type------------------
# Naive
data2 <- set.K(data,K=NULL)
params <- set.params(fixed=NULL, fixed.type=NULL)

# Populations structure by kinship and PCs
##data2 <- set.K(data)
##params <- set.params(n.PC=5, fixed=NULL, fixed.type=NULL)

#params <- set.params(fixed=c("Grp1","Grp2","Grp3","Grp4"), fixed.type=rep("numeric",4))

# GWAS execution
data3 = GWASpoly(data2, models=c("general","additive","1-dom", "2-dom"),traits=c("Gota"), params=params, n.core=4)

pdf ("out-agrosavia-gota-Naive-plots.pdf") 
	# QQ-plot Output
	message (">>>> QQ-plot...")
	par(mfrow=c(2,3)) #specifies a 2 x 3 panel
	models <- c("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
	for (i in 1:6) {
		qq.plot(data3,trait="Gota",model=models[i])
	}  

	# QTL Detection
	message (">>>> QTL Detection ...")
	data4 = set.threshold (data3, method="Bonferroni",level=0.05)
	get.QTL (data4)

	message (">>>> Manhattan plot...")
	# Manhattan plot Output
	par(mfrow=c(2,3)) #specifies a 2 x 3 panel
	models <- c("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
	for (i in 1:6) {
		manhattan.plot (data4, trait="Gota", model=models[i])
	}  
dev.off()
write.GWASpoly (data4, "Gota", "out-agrosavia-gota-Naive-scores.tbl", what="scores", "delim"="\t")


