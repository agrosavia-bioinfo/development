#!/usr/bin/Rscript
library (parallel)
args = commandArgs(trailingOnly = TRUE)


genotype  = "genotype-Enciso-gwaspoly-renamed.tbl"
phenotype = "phenotype-Enciso-LBlight-2014-renamed.tbl"
snps      = "snps-annotations.tbl"

gwasTypes = c("Naive", "Kinship", "Structure", "Kinship+Structure")
gwasTypes = c("Naive")

cmmList  = list()
for (type in gwasTypes) {
	cmm = sprintf ("Rscript gwas-polypipeline.R --pheno %s --geno %s --snps %s --model %s", 
				   phenotype, genotype, snps, type )
	cmmList = append (cmmList, cmm)
	print (cmm)
}

#print (cmmList)
mclapply (cmmList, system, mc.cores=4)


