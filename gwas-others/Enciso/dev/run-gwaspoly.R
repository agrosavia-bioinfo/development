#!/usr/bin/Rscript
library (parallel)
args = commandArgs(trailingOnly = TRUE)

phenotypesFilenames = c("phenotype-Enciso-LBlight-2010.tbl",
						"phenotype-Enciso-LBlight-2011.tbl",
						"phenotype-Enciso-LBlight-2012.tbl",
						"phenotype-Enciso-LBlight-2013.tbl",
						"phenotype-Enciso-LBlight-2014.tbl",
						"phenotype-Enciso-LBlight-2015.tbl",
						"phenotype-Enciso-LBlight-2017.tbl")

genotype = "genotype-Enciso-gwaspoly.tbl"
cmmList  = list()
for (phenotype in phenotypesFilenames) {
	cmm = sprintf ("Rscript gwaspoly-structure-simple.R %s %s", phenotype, genotype)
	cmmList = append (cmmList, cmm)
	print (cmm)
}

mclapply (cmmList, system, mc.cores=4)


