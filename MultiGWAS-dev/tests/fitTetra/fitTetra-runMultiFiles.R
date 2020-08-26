#!/usr/bin/Rscript
USAGE="
Run fitTetra in parallel on allele signal ratios from tetraploids samples
fitTetra-parallel.R <file of tetraploid ratios>
"
source ("lglib02.R")

library (fitTetra)
library (parallel)

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

#---------------------------------------------------------------------
# Parallel fitTetra
#---------------------------------------------------------------------
fitTetraFile <- function (inDir, outDir, fileFT) {
	message (fileFT,"..." )
	x=1
	inFile = paste0 (inDir, "/", fileFT)
	fitTetraData  = read.csv (inFile)

	df.tetra <- with(fitTetraData, data.frame(MarkerName=MarkerName, 
                     SampleName=SampleName, ratio=X_Raw/(X_Raw+Y_Raw)))

	markerName = df.tetra$MarkerName[1]
	outFile = paste0 (outDir, "/", fileFT)
	msg=""

	out = tryCatch ({
		unmix <- fitTetra(marker=x, data=df.tetra)
		samples    = unmix$scores$sample
		genos      = unmix$scores$geno 
		names (genos) = samples
		genosDF    = data.frame (Marker=markerName, t(genos))

		message ("...",fileFT)
		write.csv (genosDF, outFile, quote=F, row.names=F)
		msg = paste0 (fileFT, ",OK")
		return (msg)
	},
	error=function (cond) {
		message ("Error in ", fileFT)
		msg = paste0 (fileFT, ",Error")
		return (msg)
	})
	return (out)
}
#---------------------------------------------------------------------
#---------------------------------------------------------------------

#fitTetraData = read.csv ("CCC_Andigena_677_2015_fitTetra.csv")
inDir = "in-splits"
outDir = "out-fitTetra"
createDir (outDir)

files = list.files (inDir)

#fitTetraData  = read.csv (inFile)
#mx = fitTetraX (2, df.tetra)

fitTetraFile (inDir, outDir, files[[2]])

out=mclapply (files, fitTetraFile, inDir=inDir, outDir=outDir, mc.cores=8)
write.table (unlist (out), "log-fitTetra.csv", quote=F, row.names=F, col.names=F, sep=",")

