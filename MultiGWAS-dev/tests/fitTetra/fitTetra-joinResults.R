#!/usr/bin/Rscript
USAGE="
Join the results from fitTetra by SNP
INPUT: directory of individual results
OUTPUT: a single results table of genotypes"

inputDir = "out-fitTetra"

inputFiles = list.files (inputDir)

resultsTable = data.frame ()

for (f in inputFiles) {
	message (f, "...")
	inFile = paste0 (inputDir, "/", f)
	fileTable = read.csv (inFile)
	resultsTable = rbind (resultsTable, fileTable)
}
write.csv (resultsTable, "out-results-fitTetra.csv", quote=F, row.names=F)
