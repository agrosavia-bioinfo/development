#!/usr/bin/python
#------------------------------------------------------------------
# Pipeline to run GWAS using pling for a input genotype and phenotype
# r0.1: Python modules with many confused input parameters
#       It works only for the original example
#------------------------------------------------------------------
import sys, os

args = sys.argv

#args = ["","atwell2010_geno-20.ped", "atwell2010_geno.map", "atwell2010_ftfield.txt", "out"]
args = ["","atwell2010_geno.ped", "atwell2010_geno.map", "atwell2010_ftfield.txt", "out"]

GENOTYPE  = args [1]
SNPSMAP   = args [2]
PHENOTYPE = args [3]
OUTDIR    = args [4]

MAXPCS    = 3
#MAXPCS    = 3

# Script to run recoding from PLINK plain text to binary file format

# Set the directory of the data
##data_dir=/home/gwasuser/data/a.thaliana

# Set the directory of the output (which in this case will correspond to input directory)
##output_dir=/home/gwasuser/data/a.thaliana

# Call to PLINK to transform atwell2010_geno_new.[ped,map] to atwell2010_geno_new.[bed,bim,fam] with the following  flags:
# --file: full path prefix to plain PLINK files, e.g. when {filename} is given, PLINK searches for {filename}.ped and {filename}.map
# --make-bed: invokes generation of the output in binary PLINK file format
# --out: full path prefix to output filename, e.g. when {filename} is given, PLINK will store the output in {filename}.bed, {filename}.bim  and {filename}.fam

#----------------------------------------------------------
#----------------------------------------------------------
def main ():

    (inGenotype, binaryGenotype) = encodeGenotypeToBinaryFormat (GENOTYPE)
    print "\n>>>", inGenotype, binaryGenotype

    filteredGenotype = preprocessGenotype (binaryGenotype)
    print "\n>>>", filteredGenotype

    # Population Structure Correction
    outEigensoft = calculatePCs (filteredGenotype, MAXPCS)
    print "\n>>>", outEigensoft

    outFileEigenvec = createCovariatesFromPCs (filteredGenotype, outEigensoft)
    print "\n>>>", outFileEigenvec

"""
    smallestInflation, outPCsPhenotype =  searchLeadingNumberOfPCs (filteredGenotype, PHENOTYPE, outFileEigenvec, MAXPCS, outEigensoft)
    print "\n>>>", smallestInflation, outPCsPhenotype
    

    outGWAS = runUnivariateGWASWithPopulationCorrection (smallestInflation, filteredGenotype, PHENOTYPE, outFileEigenvec, outPCsPhenotype)
    print "\n>>>", outGWAS

    outManhattan, outQQ = plotResultsUniGWAS (outGWAS)
    print "\n>>>", outManhattan, outQQ
"""

#----------------------------------------------------------
# Plot results Univariate GWAS with population structure correction
#----------------------------------------------------------
def plotResultsUniGWAS (outGWAS):
    msg = "Plot results Univariate GWAS with population structure correction"

    # Manhattan plot
    manhattan_cmdl = "/home/gwasuser/tools/cmdl_wrappers/manhattan_cmdl"
    inPlinkResults = "%s.assoc.linear" % (outGWAS) 
    outManhattan   = "%s_%s" % (outGWAS, "manhattan.png")
    inMethod       = "plink" 

    cmdMnhPlot = "python %s --pvalue_file %s --out_file ./%s --method %s " % \
                 (manhattan_cmdl, inPlinkResults, outManhattan, inMethod)

    execCommand (cmdMnhPlot, "Manhattan plot...")

    # QQ plot
    qq_cmdl    = "/home/gwasuser/tools/cmdl_wrappers/qq_cmdl"
    outQQ      = "%s_%s" % (outGWAS, "qq.png")
    inMethod   = "plink" 
    inLogfile  = "%s.log" % (outGWAS)

    cmdQQPlot  = "python %s --pvalue_file %s --out_file ./%s --method %s --log_file %s " % \
                 (qq_cmdl, inPlinkResults, outQQ, inMethod, inLogfile)
                 
    execCommand (cmdQQPlot, "QQ plot...")

    return outManhattan, outQQ

#----------------------------------------------------------
# Univariate GWAS with PCs to correct for population structure"
#----------------------------------------------------------
def runUnivariateGWASWithPopulationCorrection (smallestInflation, filteredGenotype, PHENOTYPE, outFileEigenvec, outPCsPhenotype):
    inGenotype   = filteredGenotype
    inCovariates = outFileEigenvec
    inPhenotype  = PHENOTYPE
    outGWAS      = "%s_%s_pcs" % (outPCsPhenotype, smallestInflation)

    msg = "Univariate GWAS with PCs to correct for population structure"
    cmdStrPC2 = "plink --bfile %s --out %s --pheno %s --linear --covar %s --covar-number 1-%s --allow-no-sex --adjust --hide-covar" 
    cmd = cmdStrPC2 % (inGenotype, outGWAS, inPhenotype, inCovariates, smallestInflation)
    execCommand (cmd, msg)

    return outGWAS


#----------------------------------------------------------
# Search for the leading number of PCs to include as covariates 
# for Univariate GWAS with PCs to correct for population structure
#----------------------------------------------------------
def searchLeadingNumberOfPCs (filteredGenotype, PHENOTYPE, outFileEigenvec, MAXPCS, outEigensoft):
    cmd = "Script to correct for population structure using principal components (PCs) as covariates in PLINK"
    inGenotype   = filteredGenotype
    inPhenotype  = PHENOTYPE
    inCovariates = outFileEigenvec

    namePhenotype   = PHENOTYPE.split (".")[0].split("_")[1]
    outPCsPhenotype = outEigensoft + "_" + namePhenotype

    mp = 0
    cmdStrPC0 = "plink --bfile %s --out %s_%s_pcs --pheno %s --linear --allow-no-sex --adjust --hide-covar" 
    cmdStrPC1 = "plink --bfile %s --out %s_%s_pcs --pheno %s --linear --covar %s --covar-number 1 --allow-no-sex --adjust --hide-covar" 
    cmdStrPC2 = "plink --bfile %s --out %s_%s_pcs --pheno %s --linear --covar %s --covar-number 1-%s --allow-no-sex --adjust --hide-covar" 

    logFilesList = []
    for mp in range (MAXPCS):
        if   (mp == 0): cmd = cmdStrPC0 % (inGenotype, outPCsPhenotype, mp, inPhenotype)
        elif (mp == 1): cmd = cmdStrPC1 % (inGenotype, outPCsPhenotype, mp, inPhenotype, inCovariates)
        else:           cmd = cmdStrPC2 % (inGenotype, outPCsPhenotype, mp, inPhenotype, inCovariates, mp)

        print ">>> plink with ", mp, "PCs"
        execCommand (cmd)

        cmds = ["rm %s_%s_pcs.%s" % (outPCsPhenotype, mp, x) for x in ["assoc.linear","assoc.linear.adjusted","nosex"]]
        for i in cmds: execCommand (i)

        logFile = "%s_%s_pcs.log" % (outPCsPhenotype, mp) 
        logFilesList.append (logFile)

    print ">>>>>>>>>>>>>>>>>>>>"
    for i in logFilesList:
        print i

    smallestInflation = searchLowValueFromPCsLogs (logFilesList)

    return smallestInflation, outPCsPhenotype

#----------------------------------------------------------
# Parse de PCs file to create covariates
#----------------------------------------------------------
def createCovariatesFromPCs (filteredGenotype, outEigensoft):
    msg = "parse the file such that it can be used as a covariate file for PLINK or FaST-LMM"
    parser_cmdl="/home/gwasuser/tools/cmdl_wrappers/parse_eigensoft_cmdl"
    inFileEigenvec  = outEigensoft + ".eigenvecs.evec"
    inFileFam       = filteredGenotype + ".fam"
    outFileEigenvec = outEigensoft + ".eigenvec"
    cmd = "python %s --eigenvec_file %s --fam %s --out %s" % \
            (parser_cmdl,  inFileEigenvec, inFileFam, outFileEigenvec)
    execCommand (cmd, msg)

    msd = "Delete unnecessary files..."
    cmd = " rm %s.eigenvals %s.eigenvecs.par %s.eigenvecs" % (outEigensoft, outEigensoft, outEigensoft)
    execCommand (cmd, msg)

    return outFileEigenvec


#----------------------------------------------------------
# Calculate PCAs for Population Structure Correction
#----------------------------------------------------------
def calculatePCs (filteredGenotype, MAXPCS):
    msg = "Script to obtain the covariates -- first %s principal components (PCs) using EIGENSOFT" % MAXPCS
    inFile  = filteredGenotype
    outEigensoft = filteredGenotype + "_eigensoft"
    cmm = "smartpca.perl -i %s.bed -a %s.bim -b %s.fam -k %s -o %s.eigenvecs -p %s.pca_plot -l %s.log -e %s.eigenvals " % \
           (inFile, inFile, inFile, MAXPCS, outEigensoft, outEigensoft, outEigensoft, outEigensoft)
    execCommand (cmm, msg)
    return (outEigensoft)

#----------------------------------------------------------
# Preprocess by filtering the genotype
#----------------------------------------------------------
def preprocessGenotype (binaryGenotype):
    msg = "Script to run data preprocessing with PLINK"
    filteredGenotype = binaryGenotype + "_filtered"
    cmm = "plink --bfile %s --out %s --make-bed --hwe 1e-100 --maf 0.01 --mind 0.1 --geno 0.1" % (binaryGenotype, filteredGenotype)
    execCommand (cmm, msg)
    return (filteredGenotype)

#----------------------------------------------------------
# Enconding files from text to binary format
#----------------------------------------------------------
def encodeGenotypeToBinaryFormat (genotype):
    msg = "Script to run recoding from PLINK plain text to binary file format"
    inGenotype  = GENOTYPE.split (".")[0] 
    binaryGenotype = inGenotype + "_bin"
    cmm = "plink --file %s --make-bed --out %s" % (inGenotype, binaryGenotype)
    execCommand (cmm, msg)
    return (inGenotype, binaryGenotype)

#----------------------------------------------------------
# Search for the lowest value from log files generated in PCAs
#----------------------------------------------------------
def searchLowValueFromPCsLogs (logFilesList):
    inflationList = []
    for logFile in logFilesList:
        lines = open (logFile).readlines ()
        for line in lines:
            if "inflation" in line:
                inflation = line.split("=")[1].strip()
                inflationList.append (inflation)

    #print (inflationList)
    smallest = inflationList.index (min (inflationList))
    return smallest

#----------------------------------------------------------
# Execute a command and print the commands
#----------------------------------------------------------
def execCommand (cmm, msg=None):
    if (msg==None):
        print (cmm)
        os.system (cmm)
    else:
        print ("-----------------------------------------------------")
        print ">>>", msg
        print ("-----------------------------------------------------")
        print (cmm)
        os.system (cmm)




#----------------------------------------------------------
#----------------------------------------------------------
main ()


