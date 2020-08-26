#!/usr/bin/Rscript
USAGE="
Run fitTetra in parallel on allele signal ratios from tetraploids samples
fitTetra-parallel.R <file of tetraploid ratios>
"
source ("lglib02.R")

library (fitPoly)
library (parallel)
USAGE="
Run fitpoly on agrosavia clustercall data
INPUT: file with polyploids samples; MarkerName, SampleName, and ratio
OUTPUT: A matrix with Markers vs. Samples"

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#fitData  = read.csv ("datasets/data-fitPoly-tetra.csv")
outPrefix = "out-fitPoly"

fitData = read.csv ("datasets/CCC_Andigena_677_2015_fitTetra.csv")
df.tetra <- with(fitData, data.frame(MarkerName=MarkerName, 
                 SampleName=SampleName, ratio=X_Raw/(X_Raw+Y_Raw)))
#fp=fitOneMarker (ploidy=4, marker="mrk002", data=fitData);fp
saveMarkerModels (ploidy=4, data=df.tetra, filePrefix=outPrefix,
				  allModelsFile=T, ncores=7)
