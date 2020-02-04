#!/usr/bin/Rscript
# r0.1: Functions + write to RData cache

library (GWASpoly)
args = commandArgs(trailingOnly = TRUE)

#-------------------------------------------------------------
#args = c (phenotypeFile = "agrosavia-phenotype-gwaspoly.tbl", genotypeFile  = "agrosavia-genotype-gwaspoly-checked.tbl")
data = data2 = data3 = NULL
testModels = snpModels = testTraits = nTraits = gwasModel = gwasParams = NULL
phenotype = genotype = structure = NULL

phenoFile  = "phenotype-checked.tbl"
genoFile   = "genotype-checked.tbl"
structFile = "structure-checked.tbl"
  
main <- function (args) {
	args = c ("agrosavia-phenotype-gwaspoly.tbl", 
			  "agrosavia-genotype-gwaspoly-checked.tbl",
			  "K5_estructura_tetraploides_2017.txt")

	print (args)
	phenotypeFile = args [1]
	genotypeFile  = args [2]
	structureFile = args [3]

	# Check matchs in files (samples, samples) and write new files
	data = checkFiles (phenotypeFile, genotypeFile, structureFile)
	phenotype <<- data [[1]]
	genotype  <<- data [[2]]
	structure <<- data [[3]]

	# Load data
	if (file.exists ("gwas.RData")) load (file="gwas.RData")
	msg (ls())

	# Set initials for the gwas: traits, models, nTraits
	initials = setInitials () 
	testModels <<- initials [[1]]
	snpModels <<- initials [[2]]
	testTraits <<- initials [[3]]
	nTraits    <<- initials [[4]]

	# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
	data1 <<- readData (phenoFile, genoFile, nTraits, data1)

	gwasModel = "Naive"

	# Set the kinship
	data2 <<- setKinship (data1, data2, gwasModel)

	# Populations structure and kinship
	gwasParams = setPopulationStructure (genotype, structure, gwasModel)

	# GWAS execution
	data3 = runGwaspoly (data2, snpModels, testTraits, gwasParams, data3)

	# Plot results
	plotResults (data3, testModels, testTraits, gwasModel)

	save(data, data2, data3, file="gwas.RData") 
}

#-------------------------------------------------------------
#-------------------------------------------------------------
checkFiles <- function (phenotypeFile, genotypeFile, structureFile) {
	msg ("Checking files...")
	phenotypeAll = read.table (phenotypeFile, header=T, sep=",")
	genotypeAll  = read.table (genotypeFile, header=T, sep=",")
	structureAll = read.table (structureFile, header=T)

	samplesPheno  = phenotypeAll$sample
	samplesGeno   = colnames (genotypeAll)
	samplesStruct = structureAll$sample

	commonMarkers = Reduce (intersect, list (samplesGeno, samplesPheno, samplesStruct))

	phenotype  = phenotypeAll [phenotypeAll$sample %in% commonMarkers,]
	genoColums = c("Markers","Chrom","Position", commonMarkers)
	genotype   = genotypeAll  [,colnames(genotypeAll) %in% genoColums]
	structure  = structureAll [structureAll[,1] %in% commonMarkers,]

	write.table (file="phenotype-checked.tbl", phenotype, row.names=F, quote=F, sep=",")
	write.table (file="genotype-checked.tbl", genotype, row.names=F, quote=F, sep=",")
	write.table (file="structure-checked.tbl", structure, row.names=F, quote=F, sep=",")

	return (list (phenotype, genotype, structure))

}
#-------------------------------------------------------------
# Define a set of initial models, traits, and gwas analyses
#-------------------------------------------------------------
setInitials <- function () {
	msg ("Setting initials...")
	#testTraits  = c ("tuber_shape")
	#testModels <- c ("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
	#testTraits  = c ("tuber_shape")
	#testModels  = snpModels = c ("general")
	testModels  = snpModels = c ("additive", "general")
	testTraits <- c ("gota")
	msg (">>>> Models: ", testModels)
	msg (">>>> Traits: ", testTraits)
	nTraits=1

	return (list (testModels, snpModels, testTraits, nTraits))
}

#-------------------------------------------------------------
# Set the kinship matrix using the realized relationship matrix  
#-------------------------------------------------------------
setKinship <- function (data, data2, gwasModel) {
	if (!is.null (data2)) {  # When data is previously loaded
		msg ("Loading kinship...")
		return (data2)
	}

	kinshipMatrix = NULL

	if (gwasModel == "Naive") {
		msg ("\t>>>>", "Naive model: noKnoQ") 
		markerNames     = data@pheno [,1]
		n               = length (markerNames)
		kinshipMatrix   = matrix (diag (n), n, n, dimnames=list (markerNames, markerNames))
		data2           = set.K (data, K=kinshipMatrix)
	}else {
		msg ("Calculating kinship...")
		data2           = set.K(data)
	}		
	return (data2)
}
		
#-------------------------------------------------------------
# Fix Populations structure and kinship
#-------------------------------------------------------------
setPopulationStructure <- function (genotype, structure, gwasModel) {
	if (gwasModel=="Structure") {
		#genotypeStructure = cbind 
	}


	#params <- set.params(fixed=c("Grp1","Grp2","Grp3","Grp4"), fixed.type=rep("numeric",4))
	#params <- set.params(n.PC=10)
	gwasParams = NULL

	return (gwasParams)
}

#-------------------------------------------------------------
# GWAS execution
#-------------------------------------------------------------
runGwaspoly <- function (data2, snpModels, testTraits, gwasParams, data3) {
	if (!is.null (data3)) {  # When data is previously loaded
		msg ("Loading GWASpoly...")
		return (data3)
	}

	msg ("Running GWASpoly...")

	#snpModels = c("general","additive","1-dom", "2-dom")
	data3 = GWASpoly(data2, models=snpModels, traits=testTraits, params=gwasParams)

	return (data3)
}

#-------------------------------------------------------------
# Plot results
#-------------------------------------------------------------
plotResults <- function (data3, testModels, testTraits, gwasModel) {
	msg (">>> Plotting results...")

	# QTL Detection
	data4 = set.threshold (data3, method="Bonferroni",level=0.05)
	get.QTL (data4)

	# QQ-plot Output
	  #par(mfrow=c(2,3)) #specifies a 2 x 3 panel
	for (i in 1:length(testModels)) {
		plotName = sprintf("plot-gwaspoly-qq-%s-%s.pdf", gwasModel, testModels[i] )
		pdf (file=plotName, width=7,height=7)
		qq.plot(data3,trait=testTraits [1], model=testModels[i])
		dev.off()

		plotName = sprintf("plot-gwaspoly-manhattan-%s-%s.pdf", gwasModel, testModels[i] )
		pdf (file=plotName, width=7,height=7)
		manhattan.plot (y.max=20,data4, trait=testTraits[1], model=testModels [i])
		dev.off()
	}

	# Manhattan plot Output
	# pdf (file="plots-gwaspoly-manhattan.pdf")
	 #  #par(mfrow=c(2,3)) #specifies a 1 x 3 panel
	  # for (i in 1:length (gwasModels)) {
	# 	  text ("Naive")
	# 	  for (j in 1:length(testModels)) {
	# 		manhattan.plot (y.max=20,data4, trait=testTraits[1], model=testModels [j])
	# 	  }  
	 #  }
	#dev.off ()
}
#-------------------------------------------------------------
# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
#-------------------------------------------------------------
readData <- function (phenotypeFile, genotypeFile, nTraits, data1) {
	msg ("Reading data...")

	if (!is.null (data))  # When data is previously loaded
		return (data)

	data1 = read.GWASpoly (ploidy = 4, pheno.file = phenotypeFile, geno.file = genotypeFile, 
						  format = "numeric", n.traits = nTraits, delim=",")

	return (data1)
}

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) {
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}

#-------------------------------------------------------------
# Call main 
#-------------------------------------------------------------
main ()


