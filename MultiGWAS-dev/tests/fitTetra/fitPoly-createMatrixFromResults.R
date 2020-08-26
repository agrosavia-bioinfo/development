#!/usr/bin/Rscript
USAGE="
Create a matrix of Markers vs Samples from fitpoly results
INPUT: Filtered fitPoly results: MarkerName, SampleName, geno
       start and end markers (integers), label of output 
OUTPUT: File with a Matrix of Markers vs Samples
"

source ("lglib02.R")
library (parallel)

createMatrix <- function (snp, fitGenos) {
	message ("...", snp)
	dataSnp = fitGenos [fitGenos$MarkerName==snp,]
	genos = data.frame (t(dataSnp$geno))
	if (all (is.na (genos))) 
		return (NULL)
	names (genos) = dataSnp$SampleName
	rownames (genos) = snp
	return (genos)
 }

#fitResultsFile = "results/out-fitPoly_scores.dat"
#fitResults = read.table (fitResultsFile, header=T, sep="\t")

#fitGenos     = fitResults [,c("MarkerName","SampleName","geno")]
fitGenosFile = "results/out-fitPoly_scores-genos.dat"
#write.csv (fitGenos, fitGenosFile, quote=F, row.names=F)

#----
args = commandArgs(trailingOnly = TRUE)
start = args [1]
end   = args [2]
label = args [3]

fitGenos  = read.csv (fitGenosFile)		
snpList    = levels (fitGenos$MarkerName)

outs = mclapply (snpList[start:end], createMatrix, fitGenos, mc.cores=7)
fitGenosDF = do.call (rbind.data.frame, outs)

outFilename = addLabel ("out-fitPoly-matrix.csv", label)

fitPolyMatrix = data.frame (Makers=rownames (fitGenosDF), fitGenosDF)
write.csv (fitPolyMatrix, outFilename, quote=F, row.names=F)


