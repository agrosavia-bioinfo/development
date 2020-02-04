#!/usr/bin/Rscript
library (GWASpoly)
args = commandArgs(trailingOnly = TRUE)

#-------------------------------------------------------------
#args = c (phenotypeFile = "agrosavia-phenotype-gwaspoly.tbl", genotypeFile  = "agrosavia-genotype-gwaspoly-checked.tbl")
  
main <- function (args) {
	data = data3 = popResults = NULL
	args = c ("agrosavia-phenotype-gwaspoly.tbl", "agrosavia-genotype-gwaspoly-checked.tbl")
	print (args)
	phenotypeFile = args [1]
	genotypeFile  = args [2]

	# Load data
	if (file.exists ("gwas.RData")) load (file="gwas.RData")
	msg (ls())

	# Set initials for the gwas: traits, models, nTraits
	initials = setInitials () 
	testModels = initials [[1]]
	gwasModels = initials [[2]]
	testTraits = initials [[3]]
	nTraits    = initials [[4]]

	# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
	data = readData (phenotypeFile, genotypeFile, nTraits, data)

	# Populations structure and kinship
	popResults = fixPopulationStructure (data, popResults)
	data2      = popResults [[1]]
	gwasParams = popResults [[2]]

	# GWAS execution
	data3 = runGwaspoly (data2, gwasModels, testTraits, gwasParams, data3)

	# Plot results
	plotResults (data3, testModels, testTraits)

	save(data, data3, popResults, file="gwas.RData") 
}

#-------------------------------------------------------------
#-------------------------------------------------------------
setInitials <- function () {
	msg ("Setting initials...")
	#testTraits  = c ("tuber_shape")
	#testModels <- c ("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
	#testTraits  = c ("tuber_shape")
	testModels  = gwasModels = c ("additive")
	testTraits <- c ("gota")
	msg (">>>> Models: ", testModels)
	msg (">>>> Traits: ", testTraits)
	nTraits=1

	return (list (testModels, gwasModels, testTraits, nTraits))
}

#-------------------------------------------------------------
# Fix Populations structure and kinship
#-------------------------------------------------------------
fixPopulationStructure <- function (data, popResults) {
	if (!is.null (popResults)) {  # When data is previously loaded
		msg ("Loading kinship...")
		data2      = popResults [[1]]
		gwasParams = popResults [[2]]
		return (list (data2, gwasParams))
	}
	msg ("Calculating kinship...")
	data2 <- set.K(data)

	#params <- set.params(fixed=c("Grp1","Grp2","Grp3","Grp4"), fixed.type=rep("numeric",4))
	gwasParams = NULL
	#params <- set.params(n.PC=10)

	return (list (data2, gwasParams))
}

#-------------------------------------------------------------
# GWAS execution
#-------------------------------------------------------------
runGwaspoly <- function (data2, gwasModels, testTraits, gwasParams, data3) {
	if (!is.null (data3)) {  # When data is previously loaded
		msg ("Loading GWASpoly...")
		return (data3)
	}

	msg ("Running GWASpoly...")

	#gwasModels = c("general","additive","1-dom", "2-dom")
	data3 = GWASpoly(data2, models=gwasModels, traits=testTraits, params=gwasParams)

	return (data3)
}

#-------------------------------------------------------------
# Plot results
#-------------------------------------------------------------
plotResults <- function (data3, testModels, testTraits) {
	msg (">>> Plotting results...")
	# QQ-plot Output
	pdf (file="plots-qq-gwaspoly.pdf")
	  par(mfrow=c(2,3)) #specifies a 2 x 3 panel
	  models <- testModels
	  for (i in 1:length(testModels)) {
		qq.plot(data3,trait=testTraits [1], model=models[i])
	  }
	dev.off()

	# QTL Detection
	data4 = set.threshold (data3, method="Bonferroni",level=0.05)
	get.QTL (data4)

	# Manhattan plot Output
	pdf (file="plots-manhattan-gwaspoly.pdf")
	  par(mfrow=c(2,3)) #specifies a 1 x 3 panel
	  models <- testModels
	  for (i in 1:length(testModels)) {
		manhattan.plot (data4, trait=testTraits[1], model=models [i])
	  }  
	dev.off ()
}
#-------------------------------------------------------------
# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
#-------------------------------------------------------------
readData <- function (phenotypeFile, genotypeFile, nTraits, data) {
	msg ("Reading data...")

	if (!is.null (data))  # When data is previously loaded
		return (data)

	data = read.GWASpoly (ploidy = 4, pheno.file = phenotypeFile, geno.file = genotypeFile, 
						  format = "numeric", n.traits = nTraits, delim=",")

	return (data)
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


