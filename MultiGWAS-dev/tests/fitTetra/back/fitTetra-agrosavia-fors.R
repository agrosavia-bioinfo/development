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
# Single marker, multiple mixture models
fitTetraX <- function (x, df.tetra) {
	print (as.character (df.tetra$MarkerName[x]))
	unmix <- fitTetra(marker=x, data=df.tetra)
	#if (is.na (unmix$scores)) return (NA)

	SNP    = as.character (df.tetra$MarkerName [x])
	genos  = unmix$scores$geno 
	rowg   = c (SNP, as.character (genos))
	return (rowg)
}

fitTetraFile <- function (inDir, outDir, fileFT) {
	message (fileFT,"..." )
	x=1
	inFile = paste0 (inDir, "/", fileFT)
	fitTetraData  = read.csv (inFile)

	df.tetra <- with(fitTetraData, data.frame(MarkerName=MarkerName, 
                     SampleName=SampleName, ratio=X_Raw/(X_Raw+Y_Raw)))

	unmix <- fitTetra(marker=x, data=df.tetra)
	if (is.na (unmix$scores)) 
		return (NA)

	markerName = df.tetra$MarkerName[1]
	samples    = unmix$scores$sample
	genos      = unmix$scores$geno 
	names (genos) = samples
	genosDF    = data.frame (Marker=markerName, t(genos))

	outFile = paste0 (outDir, "/", fileFT)
	message ("...",fileFT)
	write.csv (genosDF, outFile, quote=F, row.names=F)
}
#---------------------------------------------------------------------
#---------------------------------------------------------------------

#fitTetraData = read.csv ("CCC_Andigena_677_2015_fitTetra.csv")
inDir = "out-splits"
outDir = "out-fitTetra"
createDir (outDir)

files = list.files (inDir)

#fitTetraData  = read.csv (inFile)
#mx = fitTetraX (2, df.tetra)

#fitTetraFile (inDir, outDir, files[[3]])

mclapply (files, fitTetraFile, inDir=inDir, outDir=outDir, mc.cores=1)

for (f in files) {
	message ("...", f)
	inFile = paste0 (inDir, "/", f)
	fitTetraData  = read.csv (inFile)
	colnames (fitTetraData) = c("MarkerName", "SampleName", "X_Raw","Y_Raw","Theta","R")

	df.tetra <- with(fitTetraData, data.frame(MarkerName=MarkerName, 
                     SampleName=SampleName, ratio=X_Raw/(X_Raw+Y_Raw)))

	N = length (df.tetra$MarkerName)

	allList  = mclapply (1:N, FUN=fitTetraX, df.tetra, mc.cores=1)
	allTrue  = allList [which (!is.na (allList))]

	if (length (allTrue) > 0) {
		allDF =data.frame (matrix (unlist (allTrue), nrow=length (allTrue), byrow=T))
		sampleNames = as.character (unique (df.tetra$SampleName))
		colnames (allDF) = c("Marker", sampleNames)
		outFile = paste0 (outDir, "/", f)
		write.csv (allDF, outFile, quote=F, row.names=F)
	}else {
		message (paste0 ("...Empty ", f))
		logFile = file ("empty-files.txt", open="a")
		writeLines (f, logFile)
		close (logFile)
	}
}

