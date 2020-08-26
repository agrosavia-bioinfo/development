#!/usr/bin/Rscript
source ("lglib02.R")

library (fitTetra)
library (parallel)

#--------------------------------------------------------
# Split the files of an input directory in bins according to the\n
# size of the bin. The bins are put in an output directory\n
#
# outputDir is the destiny dir for bins
# binSize is the number of file by bin
# sizeFill is the prefix for each bin filename
#--------------------------------------------------------
createBins <- function (inputData, outputDir) {
	createDir (outputDir)

	dataFT  = read.csv (inputData)
	markers = as.character (unique (dataFT$MarkerName))

	for (k in markers) {
		message ("...", k)
		dataFTM = dataFT [dataFT$MarkerName %in% k,]
		binFilename = paste0 (outputDir, "/SNP_", k, ".csv")
		write.csv (dataFTM, binFilename, quote=F, row.names=F)
	}
}

#----------------------------------------------------------
# Create dir, if it exists the it is renamed old-XXX
#----------------------------------------------------------
createDir <- function (newDir) {
	checkOldDir <- function (newDir) {
		name  = basename (newDir)
		path  = dirname  (newDir)
		if (dir.exists (newDir) == T) {
			oldDir = sprintf ("%s/old-%s", path, name)
			if (dir.exists (oldDir) == T) {
				checkOldDir (oldDir)
			}

			file.rename (newDir, oldDir)
		}
	}

	checkOldDir (newDir)
	system (sprintf ("mkdir %s", newDir))
}

#----------------------------------------------------------
#----------------------------------------------------------


fitTetraFile = "CCC_Andigena_677_2015_fitTetra.csv"
#fitTetraData = read.csv ("test100.csv")

createBins (fitTetraFile, "out-splits")
