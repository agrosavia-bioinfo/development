#!/usr/bin/Rscript
# r2.5: Generating two tables: bestN and significatives
# r2.41: Changed to parameters by config file. Beginning changes...
# r2.3: Removed annotations as it fails. Before config by file instead parameters
# r2.2: Write model type and inflation factor to significative QTLs

N_QTLs = 20     # Max. Number of best significative QTLs to report
args = commandArgs(trailingOnly = TRUE)
args = c("naive-4g-ACGT.config")

USAGE="USAGE: Rscript gwas-polypiline.R <config file>"
if (length (args) != 1) {
	message (USAGE)
	quit()
}

library (GWASpoly)                  # For GWAS
suppressMessages (library (dplyr))  # For data handing
suppressMessages(library (config))  # For read config file

setClass ("GWASpolyStruct", slots=c(params="list"),#testModels="character" snpModels="character", 
		  contains="GWASpoly.K")

options (width=300)

#gwasModelsTypes        = c("Naive", "Kinship", "Structure", "PCs", "Kinship+Structure", "Kinship+PCs", "Structure+PCs", "Kinship+Structure+PCs")
#snpModels              = c("general","additive","1-dom", "2-dom")
multipleTestingModel    = "Bonferroni" # Bonferroni | FDR | permute
#snpModels = testModels = c("general")
#-------------------------------------------------------------
# Global configs
#-------------------------------------------------------------
#genotypeFormati = "numeric"|"AB"|"ACGT"
#genotypeFormat = "numeric"
#gwasModel = "Naive"|"Structure"|"Kinship"|"Kinship+Structure|PCs"
#gwasModel      = "Kinship"
#PLOIDY         = 4
#-------------------------------------------------------------
data = data1 = data2 = data3 = data4 = NULL
phenotype = genotype = structure = NULL

USAGE="USAGE: Rscript gwas-polypiline.R <config file>"
main <- function (args) 
{
	params = config::get (file=args [1])

	# Read and check command line arguments
	params = setGwaspolyModels (params)

	msg (">>>>>>>>>>>>", params$gwasModel, "<<<<<<<<<<<") 

	# Load cache data
	if (file.exists ("gwas.RData")) load (file="gwas.RData")

	# Read and check matchs in params (phenotype and genotype) and write new params
	data = readCheckFiles (params$phenotypeFile, params$genotypeFile, params$structureFile)

	# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
	data1 = initGWAS (params$phenotypeFile, params$genotypeFile,
					   params$genotypePloidy, params$genotypeFormat, data1)

	# Set the kinship
	data2 <- setKinship (data1, params$gwasModel, data2)

	# Populations structure and kinship
	data3 <- setPopulationStructure (data2, params$gwasModel, data$phenotype, 
									data$structure, data3)

	# GWAS execution
	data4 <- runGwaspoly (data3, params$gwasModel, params$snpModels, data4)

	save(data, data1, data2, data3, data4, file="gwas.RData") 

	# Plot results
	if (params$genotypePloidy==4) uploidyLabel = "Tetra" else ploidyLabel = "Diplo"
	showResults (data4, params$testModels, data$trait, params$gwasModel, 
	             params$phenotypeFile, params$annotationsFile, params$genotypePloidy)

}
#-------------------------------------------------------------
# Read and check command line arguments
# Currently: phenotye, genotype, structure
#-------------------------------------------------------------
setGwaspolyModels <- function (params) 
{
	msg("Reading Files...")

	if (params$genotypePloidy=="4") {
		#testModels  = snpModels = c ("general")
		snpModels  = c("general","additive","1-dom", "2-dom")
		testModels = c("general", "additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")

		params = append (params, list (snpModels  = snpModels))
		params = append (params, list (testModels = testModels))
	}else{
		params = append (params, list (snpModels  = c("general","additive","1-dom")))
		params = append (params, list (testModels = c("general", "additive","1-dom-alt","1-dom-ref")))
	}

	print (params)
	message ("------------------------------------------------")
	message ("Summary of input parameters:")
	message ("------------------------------------------------")
	msg ("Phenotype filename : ", params$phenotypeFile) 
	msg ("Genotype filename  : ", params$genotypeFile) 
	msg ("Structure filename : ", params$structureFile) 
	msg ("SNPs filename      : ", params$annotationsFile) 
	msg ("GwAS model         : ", params$gwasModel) 
	msg ("Genotype format    : ", params$genotypeFormat) 
	msg ("Genotype ploidy    : ", params$genotypePloidy)
	msg ("Kinship            : ", params$kinship) 
	message ("------------------------------------------------")

	return (params)

}

#-------------------------------------------------------------
# Read and check files, sample samples
#-------------------------------------------------------------
readCheckFiles <- function (phenotypeFile, genotypeFile, structureFile) 
{
	msg();msg ("Reading and checking files (same samples for pheno, geno, and struct)...");msg()

	phenotypeAll = read.csv (phenotypeFile, header=T, sep=",", check.names=F)
	genotypeAll  <<- read.csv (genotypeFile, header=T, sep=",", check.names=F)

	if (is.null (structureFile)) structureAll <- NULL 
	else structureAll <- read.csv (structureFile, header=T,check.names=F)


	samplesPheno  = phenotypeAll[,1]
	samplesGeno   = colnames (genotypeAll)
	#samplesStruct = structureAll$sample

	commonMarkers = Reduce (intersect, list (samplesGeno, samplesPheno)) #, samplesStruct)

	phenotype  = phenotypeAll [phenotypeAll$Samples %in% commonMarkers,]
	genoColums = c("Markers","Chrom","Position", commonMarkers)
	genotype   = genotypeAll  [,colnames(genotypeAll) %in% genoColums]
	#structure  = structureAll [structureAll[,1] %in% commonMarkers,]

	write.table (file="phenotype-checked.tbl", phenotype, row.names=F, quote=F, sep=",")
	write.table (file="genotype-checked.tbl", genotype, row.names=F, quote=F, sep=",")
	#write.table (file="structure-checked.tbl", structure, row.names=F, quote=F, sep=",")

	trait  <- colnames (phenotype)[2]
	msg  (">>>> Evaluating trait ", trait)

	return (list (phenotype=phenotype, genotype=genotype,structure=structureAll, trait=trait))
}

#-------------------------------------------------------------
# Set the kinship matrix using the realized relationship matrix  
#-------------------------------------------------------------
setKinship <- function (data1,  gwasModel, data2) 
{
	msg();msg("Setting kinship...");msg()
 	# Load data instead calculate it
	if (!is.null (data2)) {msg ("Loading kinship..."); return (data2) }

	#kinshipMatrix = NULL
	if (gwasModel %in% c("Kinship", "Kinship+Structure", "Kinship+PCs")) {
		msg (">>>> With default kinship... ")
		kinshipMatrix = NULL
	}else { 
		msg (">>>> Without kinship...") 
		markerNames   = data1@pheno [,1]
		n             = length (markerNames)
		kinshipMatrix = matrix (diag (n), n, n, dimnames=list (markerNames, markerNames))
	}		
	data2  = set.K (data1, K=kinshipMatrix)
	return (data2)
}
		
#-------------------------------------------------------------
# Fix Populations structure and kinship
#-------------------------------------------------------------
setPopulationStructure <- function (data2, gwasModel, phenotype, structure, data3) 
{
	st <- structure
	msg();msg ("Setting population structure...");msg()
	#setClass ("GWASpolyStruct", slots=c(params="list"), contains="GWASpoly.K")
	data3 <- new ("GWASpolyStruct", data2)
	
	nPCs=0
	structNames = structTypes = NULL

	if (gwasModel %in% c("PCs","Kinship+PCs")) {
		nPCs=5
		msg (">>>> With nPCS=", nPCs, "...")
	}else if (gwasModel %in% c("Structure", "Kinship+Structure")) {
		msg (">>>> With population structure...")
		structNames <- colnames (structure[-1])
		structTypes <- rep ("numeric", length (structNames))

		#data3@pheno = phenoStruct
		data3@pheno = phenotype
		data3@fixed = st[,-1]
	}else 
		msg (">>>> Without populations structure")

	data3@params = set.params(n.PC=nPCs, fixed=structNames, fixed.type=structTypes)

	return (data3)
}

#-------------------------------------------------------------
# GWAS execution
#-------------------------------------------------------------
runGwaspoly <- function (data3, gwasModel, snpModels, data4) 
{
	msg();msg("Running GWASpoly...");msg()

	if (!is.null (data4)) { msg (">>>> Loading GWASpoly..."); return (data4) }

 	if (gwasModel %in% c("Naive","Kinship")) {
		msg (">>>> Without params")
		data4 = GWASpoly(data3, models=snpModels, traits=NULL, params=NULL, n.core=4)
	}else {
		msg (">>>> With params")
		data4 = GWASpoly(data3, models=snpModels, traits=NULL, params=data3@params)
	}
	
	return (data4)
}

#-------------------------------------------------------------
# Plot results
#-------------------------------------------------------------
showResults <- function (data4, testModels, trait, gwasModel, phenotypeFile, snpsAnnFile, ploidy) 
{
	msg();msg ("Plotting results...");msg()
	phenoName = strsplit (phenotypeFile, split=".tbl")[[1]][1]
	plotName = sprintf("out-gwasp%s-%s-plots.pdf", ploidy, gwasModel)

	pdf (file=plotName, width=11, height=7)
	n = length (testModels)
	
	# QQ-plot Output
	op <- par(mfrow = c(2,n), oma=c(0,0,3,0))
	for (i in 1:length(testModels)) {
		#par (cex.main=0.5, cex.lab=0.5, cex.axis=0.5, ann=T)
		qqPlot(data4,trait=trait, model=testModels[i], cex=0.3)
	}

	# QTL Detection
	data5 <<- set.threshold (data4, method=multipleTestingModel,level=0.05,n.permute=100,n.core=3)
	#significativeQTLs = get.QTL (data5)
	QTLs = getQTL (data5, snpsAnnFile, gwasModel, ploidy)

	msg (">>>>", "Writing results...")
	outFile = sprintf ("out-gwasp%s-%s-QTLs-bestN.tbl", ploidy, gwasModel) 
	write.table (file=outFile, QTLs$best, quote=F, sep="\t", row.names=F)
	outFile = sprintf ("out-gwasp%s-%s-QTLs-significatives.tbl", ploidy, gwasModel) 
	write.table (file=outFile, QTLs$significatives, quote=F, sep="\t", row.names=F)

	# Manhattan plot Output
	for (i in 1:length(testModels)) {
		#par (cex=1.5)
		manhattan.plot (y.max=20,data5, trait=trait, model=testModels [i])
	}
	plotTitle = sprintf ("GWAS %s-ploidy with %s for %s trait", ploidy, gwasModel, trait)  
	mtext(plotTitle, outer=T,  cex=1.5,  line=0)
	par(op)
	dev.off()
}

#-------------------------------------------------------------
# Extracts significant QTL
#-------------------------------------------------------------
getQTL <- function(data,snpsAnnFile, gwasModel, ploidyLabel, traits=NULL,models=NULL) 
{
	stopifnot(inherits(data,"GWASpoly.thresh"))
	map <<-data@map

	if (is.null(traits)) traits <- names(data@scores)
	else stopifnot(is.element(traits,names(data@scores)))

	if (is.null(models)) models <- colnames(data@scores[[1]])
	else stopifnot(is.element(models,colnames(data@scores[[1]])))

	msg ("Reading associations...")
	if (!is.null (snpsAnnFile)) snpsAnnotations <- read.csv (file=snpsAnnFile, header=T)

	n.model <- length(models)
	output <- data.frame(NULL)
	
	TRAIT <- traits [1]
	for (j in 1:n.model) {
		#ix <- which(data@scores[[TRAIT]][,models[j]] > data@threshold[TRAIT,models[j]])
		ix <- which (!is.na (data@scores[[TRAIT]][,models[j]])) 
		mdl <- data@scores[[TRAIT]][,models[j]] 
		markers <-  data.frame (SNP=data@map[ix,c("Marker")])
		if (!is.null (snpsAnnFile)) 
			snpAnn  <- merge (markers, snpsAnnotations, by.x="SNP",by.y="SNP_id", sort=F)[,c(2,7)]
		else
			snpAnn = "None"

		scores <- data@scores[[1]][,models[j]]
		datax = calculateInflationFactor (scores)

		n.ix <- length(ix)
		
		df = data.frame(Ploidy=rep (ploidyLabel, n.ix), Type=rep (gwasModel, n.ix),
						GC=rep(datax$delta,n.ix), Model=rep(models[j],n.ix),
						Score=round(data@scores[[TRAIT]][ix,models[j]],2),
					    Threshold=round(rep(data@threshold[TRAIT,models[j]],n.ix),2),
						Effect=round(data@effects[[TRAIT]][ix,models[j]],2),
						data@map[ix,])
						#snpAnn) 
						#stringsAsFactors=F,check.names=F)

		output <- rbind(output, df)
	}
	#out <-cbind (Type=gwasModel, output)
	#output <<- output [order(-output$GC,-output$Score),]
	output         = output %>% arrange (desc(Score)) %>% distinct (Marker, .keep_all=T) 
	best           = output %>% head(N_QTLs)
	significatives = output %>% filter (Score >=Threshold)
	return(list(best=best, significatives=significatives))
}
#-------------------------------------------------------------
# QQ plot
#-------------------------------------------------------------
qqPlot <- function(data,trait,model,cex=1,filename=NULL) 
{
	stopifnot(inherits(data,"GWASpoly.fitted"))
	traits <- names(data@scores)
	stopifnot(is.element(trait,traits))
	models <- colnames(data@scores[[trait]])
	stopifnot(is.element(model,models))
	scores <- data@scores[[trait]][,model]

	datax = calculateInflationFactor (scores)

	n <- length(datax$scores)
	unif.p <- -log10(ppoints(n))
	if (!is.null(filename)) {postscript(file=filename,horizontal=FALSE)}
	par(pty="s")
	plot(unif.p, datax$scores, pch=16,cex=cex,
		 xlab=expression(paste("Expected -log"[10],"(p)",sep="")),
		 ylab=expression(paste("Observed -log"[10],"(p)",sep="")),
		 main=paste(trait," (",model,") ",sep=""))

	mtext (bquote(lambda[GC] == .(datax$delta)), side=3, line=-2, cex=0.7)

	lines(c(0,max(unif.p)),c(0,max(unif.p)),lty=2)
	if (!is.null(filename)) {dev.off()}
	return(datax$delta)
}

#-------------------------------------------------------------
# Calculate the inflation factor from -log10 values
#-------------------------------------------------------------
calculateInflationFactor <- function (scores)
{
	remove <- which(is.na(scores))
	if (length(remove)>0) 
		x <- sort(scores[-remove],decreasing=TRUE)
	else 
		x <- sort(scores,decreasing=TRUE)

	pvalues = 10^-x
	chisq <- na.omit (qchisq(1-pvalues,1))
	delta  = round (median(chisq)/qchisq(0.5,1), 3)

	return (list(delta=delta, scores=x))
}

#-------------------------------------------------------------
# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
#-------------------------------------------------------------
initGWAS <- function (phenotypeFile, genotypeFile, ploidy, format, data1) 
{
	msg();msg ("Initializing GWAS...");msg()
	# When data is previously loaded
	if (!is.null (data)) {msg(">>>> Loading GWAS data..."); return (data)}

	data1 = read.GWASpoly (ploidy = ploidy, pheno.file = phenotypeFile, geno.file = genotypeFile, 
							  format = format, n.traits = 1, delim=",")

	return (data1)
}

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) 
{
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}

#-------------------------------------------------------------
# Call main 
#-------------------------------------------------------------
main (args)


