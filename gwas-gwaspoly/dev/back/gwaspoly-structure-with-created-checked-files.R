#!/usr/bin/Rscript
# r0.1: Functions + write to RData cache

library (GWASpoly)
args = commandArgs(trailingOnly = TRUE)
setClass ("GWASpolyStruct", slots=c(params="list"), contains="GWASpoly.K")

options (width=300)

#gwasModelsTypes = c("Naive", "Kinship", "Structure", "Kinship+Structure")
#snpModels = c("general","additive","1-dom", "2-dom")
#testModels <- c ("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
gwasModelsTypes = c("Structure")
snpModels       = c("general")
testModels      = c("general")
#gwasModelsTypes = c("Naive", "Kinship", "Structure", "Kinship+Structure")
#snpModels       = c("general","additive","1-dom", "2-dom")
#testModels      = c("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
testTraits      = c ("gota")
nTraits         = 1
#-------------------------------------------------------------
#args = c (phenotypeFile = "agrosavia-phenotype-gwaspoly.tbl", genotypeFile  = "agrosavia-genotype-gwaspoly-checked.tbl")
data = data2 = data3 = data4 = NULL
phenotype = genotype = structure = NULL

phenoFile       = "phenotype-checked.tbl"
genoFile        = "genotype-checked.tbl"
structFile      = "structure-checked.tbl"
phenoStructFile = "genostruct-checked.tbl"
  
main <- function (args) {
	args = c (phenoFile, genoFile, structFile, phenoStructFile)

	print (args)
	phenotypeFile   = args [1]
	genotypeFile    = args [2]
	structureFile   = args [3]
	phenoStructFile = args [4]

	# Read and chheck matchs in files (samples, samples) and write new files
	data = readData (phenotypeFile, genotypeFile, structureFile, phenoStructFile)
	phenotype   <- data [[1]]
	genotype    <- data [[2]]
	structure   <- data [[3]]
	phenoStruct <- data [[4]]

	# Load data
	if (file.exists ("gwas.RData")) load (file="gwas.RData")
	msg (ls())

	# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
	data1 <- initGWAS (phenoFile, genoFile, nTraits, data1)

	print (gwasModelsTypes)
	for (gwasModel in gwasModelsTypes) {
		msg (">>>>>>>>>>>>", gwasModel, "<<<<<<<<<<<") 
		# Set the kinship
		data2 <- setKinship (data1, gwasModel, data2)

		# Populations structure and kinship
		data3 <- setPopulationStructure (data2, gwasModel, phenotype, structure, phenoStruct, data3)

		# GWAS execution
		data4 <<- runGwaspoly (data3, gwasModel, snpModels, testTraits, data4)

		# Plot results
		showResults (data4, testModels, testTraits, gwasModel)
	}

	save(data, data1, data2, data3, data4, file="gwas.RData") 
}

#-------------------------------------------------------------
#-------------------------------------------------------------

readData <- function (phenotypeFile, genotypeFile, structureFile, phenoStructFile) {
	msg ("Reading files and checking colnames and rownames...")
	phenotype   = read.table (phenotypeFile, header=T, sep=",")
	genotype    = read.table (genotypeFile, header=T, sep=",")
	structure   = read.table (structureFile, header=T, sep=",")
	phenostruct = read.table (phenoStructFile, header=T, sep=",")

	return (list (phenotype, genotype, structure, phenostruct))

}

#-------------------------------------------------------------
#-------------------------------------------------------------
readCheckFiles <- function (phenotypeFile, genotypeFile, structureFile) {
	msg ("Reading files and checking colnames and rownames...")
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
	#testModels  = snpModels = c ("general")
	snpModels = c("general","additive","1-dom", "2-dom")
	testModels <- c ("additive","general","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
	#testModels  = snpModels = c ("additive", "general")
	testTraits <- c ("gota")
	msg (">>>> Models: ", testModels)
	msg (">>>> Traits: ", testTraits)
	nTraits=1

	return (list (testModels, snpModels, testTraits, nTraits))
}

#-------------------------------------------------------------
# Set the kinship matrix using the realized relationship matrix  
#-------------------------------------------------------------
setKinship <- function (data1,  gwasModel, data2) {

 	# Load data instead calculate it
	#if (!is.null (data2)) {msg ("Loading kinship..."); return (data2) }

	kinshipMatrix = NULL

	if (gwasModel %in% c("Naive", "Structure")) {
		msg ("Without kinship...") 
		markerNames     = data1@pheno [,1]
		n               = length (markerNames)
		kinshipMatrix   = matrix (diag (n), n, n, dimnames=list (markerNames, markerNames))
		data2           = set.K (data1, K=kinshipMatrix)
	}else 
		if (gwasModel %in% c("Kinship", "Kinship+Structure")) {
			msg ("Calculating kinship...")
			data2 = set.K(data1)
		}		

	return (data2)
}
		
#-------------------------------------------------------------
# Fix Populations structure and kinship
#-------------------------------------------------------------
setPopulationStructure <- function (data2, gwasModel, phenotype, structure, phenoStruct, data3) {
	#setClass ("GWASpolyStruct", slots=c(params="list"), contains="GWASpoly.K")
	data3 = new ("GWASpolyStruct", data2)

	if (gwasModel %in% c("Naive", "Kinship")) {
		msg ("Without population structure...")
		data3@params = list()
	}else 
		if (gwasModel %in% c("Structure", "Kinship+Structure")) {
			msg ("Fixing population structure...")
			structNames = colnames (structure[-1])
			structTypes = rep ("numeric", length (structNames))

			data3@pheno = phenoStruct
			data3@fixed = structure [,-1]

			data3@params = set.params(n.PC=20,fixed=structNames, fixed.type=structTypes)
		}

	#params <- set.params(fixed=c("Grp1","Grp2","Grp3","Grp4"), fixed.type=rep("numeric",4))
	#params <- set.params(n.PC=10)

	return (data3)
}
#-------------------------------------------------------------
# GWAS execution
#-------------------------------------------------------------
runGwaspoly <- function (data3, gwasModel, snpModels, testTraits, data4) {
	#if (!is.null (data4)) { msg ("Loading GWASpoly..."); return (data4) }

	msg ("Running GWASpoly...")

	#snpModels = c("general","additive","1-dom", "2-dom")
	if (gwasModel %in% c("Naive","Kinship")) {
		msg ("Running Naive GWAS...")
		data4 = GWASpoly(data3, models=snpModels, traits=testTraits, params=NULL)
	} else {
		msg ("Running GWAS with Fixed Population Structure...")
		data4 = GWASpoly(data3, models=snpModels, traits=testTraits, params=data3@params)
	}

	return (data4)
}

#-------------------------------------------------------------
# Plot results
#-------------------------------------------------------------
showResults <- function (data4, testModels, testTraits, gwasModel) {
	msg ("Plotting results...")

	# QTL Detection
	data5 = set.threshold (data4, method="Bonferroni",level=0.05)
	significativeQTLs = get.QTL (data5)

	msg (">>>>", "Writing results...")
	outFile = sprintf ("gwas-significative-QTLs-%s.tbl", gwasModel) 
	write.table (file=outFile, significativeQTLs, quote=F, sep="\t", row.names=F)

	print (significativeQTLs)

	# QQ-plot Output
	  #par(mfrow=c(2,3)) #specifies a 2 x 3 panel
	for (i in 1:length(testModels)) {
		plotName = sprintf("plot-gwaspoly-qq-%s-%s.pdf", gwasModel, testModels[i] )
		pdf (file=plotName, width=7,height=7)
		par (cex=1.5)
		qq.plot(data4,trait=testTraits [1], model=testModels[i])
		dev.off()

		plotName = sprintf("plot-gwaspoly-manhattan-%s-%s.pdf", gwasModel, testModels[i] )
		pdf (file=plotName, width=7,height=7)
		par (cex=1.5)
		manhattan.plot (y.max=20,data5, trait=testTraits[1], model=testModels [i])
		dev.off()

		#msg (">>>>", "Writing results...")
		#outFile = sprintf ("gwas-results-%s-scores.txt", gwasModel) 
		#write.GWASpoly (data5, "gota", outFile, what="scores", "delim"="\t")
		#outFile = sprintf ("gwas-results-%s-effects.txt", gwasModel) 
		#write.GWASpoly (data5, "gota", outFile, what="effects", "delim"="\t")
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
initGWAS <- function (phenotypeFile, genotypeFile, nTraits, data1) {
	msg ("Initializing GWAS...")
	# When data is previously loaded
	#if (!is.null (data)) {msg(">>>> Loading GWAS data..."; return (data)}

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


