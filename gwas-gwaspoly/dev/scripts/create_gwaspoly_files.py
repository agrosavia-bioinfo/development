#!/usr/bin/python

import os, sys

agrosaviaGenotypesFile  = "agrosavia-genotypes-original.tbl"
agrosaviaPhenotypesFile = "agrosavia-phenotype-original.tbl"
agrosaviaPhenotypeGwaspolyFile = "agrosavia-phenotype-gwaspoly.tbl"  

genomeAnnotationsFile   = "potato_8303SNPs_potato_dm_v4.03.gff3"
shipAnnotationsFile	 = "potato_infinium_8303_map_context_DM_v3_superscaffolds.txt"

#----------------------------------------------------------
#----------------------------------------------------------
def main ():
	createGenotypeForGwaspoly (agrosaviaGenotypesFile, shipAnnotationsFile, genomeAnnotationsFile)
	createPhenotypeForGwaspoly (agrosaviaPhenotypesFile, agrosaviaPhenotypeGwaspolyFile)

	#gotaPhenotypeDic	= createPhenotypeForGwaspoly (agrosaviaPhenotypesFile)
	#genomeAnnotationDic = createGenotypeGwaspoly (genomeAnnotationsFile)
	#shipArraySNPsDic	= getDicFromShip (shipAnnotationsFile)

	#createPEDFile (gotaPhenotypeDic, genomeAnnotationDic, shipArraySNPsDic)
	
#----------------------------------------------------------
#----------------------------------------------------------
def createPhenotypeForGwaspoly (agrosaviaPhenotypesFile, agrosaviaPhenotypeGwaspolyFile):
	dic = {}
	outFile	   = open (agrosaviaPhenotypeGwaspolyFile, "w")
	outFileErrors = open ("agrosavia-fenotipo-gota.errors", "w")
	
	for line in open(agrosaviaPhenotypesFile).readlines()[1:]:
		try:
			fields		 = line.split ("\t")
			sampleName	 = fields[0]
			sampleId	   = int (sampleName.split ("_")[1])
			gotaValue	  = fields[5] 
			dic [sampleId] = [sampleName, gotaValue]
		except:
			outFileErrors.write (line)

	outFile = open (agrosaviaPhenotypeGwaspolyFile, "w")
	outFile.write ("# SampletId, GotaValue\n")
	for key in sorted (dic.keys()):
		outFile.write ("%s\t%s\n" % (dic [key][0], dic[key][1]))
	outFile.close ()

	return dic

#----------------------------------------------------------
#----------------------------------------------------------
def createGenotypeForGwaspoly (agrosaviaGenotypesFile, shipAnnotationsFile, genomeAnnotationsFile):
	errorsFile = open ("errors-creating-gwaspoly-files.errors", "w")

	outFileErrors = open ("agrosavia-genome-annotations-repeatedIds.errors", "w")
	# Read genome annotations file 
	dicGenome = {}
	linesGenome = open (genomeAnnotationsFile).readlines()
	for line in linesGenome:
		fields		= line.split ()
		chromosome	= fields [0]
		snpId		= fields[-1].split("=")[2].strip(";")

		if dicGenome.get (snpId) != None:
			outFileErrors.write (line+"%s\n" % dicGenome.get (snpId))

		dicGenome [snpId] = chromosome

	# Read ship annotations file
	dicShip = {}
	linesShip = open (shipAnnotationsFile).readlines()[1:]
	curFile = shipAnnotationsFile
	for line in linesShip:
		try:
			print line
			fields   = line.split ()
			snpId	= fields [0]
			position = int (fields [2])
			dicShip [snpId] = position
		except:
			errorsFile.write (">>>%s:\n%s\n" % (curFile, line))

	# Write genotype file
	linesGenoAgrosavia  = open (agrosaviaGenotypesFile).readlines()
	curFile =  agrosaviaGenotypesFile
	outFile	  = open ("agrosavia-genotype-gwaspoly.tbl", "w")
	outFile.write ("# Marker, Chromosome, Position, Genotype----\n")
	for line in linesGenoAgrosavia:
		try:
			fields	 = line.split (",")
			snpId	  = fields [0]
			samplesIds = fields [1:]

			if snpId == "Markers":
				fields.insert (1, "Chrom")
				fields.insert (2, "Position")
			else:
				chrom	= dicGenome [snpId]
				position = dicShip [snpId]
				fields.insert (1, chrom)
				fields.insert (2, position)

			newLine = ",".join (map(str,fields))
			outFile.write (newLine)
		except:
			errorsFile.write (">>>%s:\n%s\n" % (curFile, line))


#----------------------------------------------------------
#----------------------------------------------------------
def getDicFromShip (shipAnnotationsFile):
	linesList = open (shipAnnotationsFile).readlines()[1:]
	dic = {}

	outFileErrors = open ("agrosavia-ship-annotations.errors", "w")
	for line in linesList:
		fields	  = line.split ()
		snpId	   = fields [0]
		position	= fields [2]
		sequence	= fields [3]
		alleles	 = getAlleles (sequence)

		repeated	  = dic.get (snpId)
		if repeated != None:
			outFileErrors.write (line)
			outFileErrors.write ("%s\n" % repeated)

		dic [snpId] = [position, alleles]

	outFile = open ("agrosavia-ship-annotations.tbl", "w")
	outFile.write ("# snpID, Position, Alleles\n")
	for key in sorted (dic.keys()):
		outFile.write ("%s\t%s\t%s\n" % (key, dic [key][0], dic [key][1]))

def getAlleles (sequence):
	allele1 = sequence.split ("[")[1].split("/")[0]
	allele2 = sequence.split ("/")[1].split("]")[0]

	return (allele1+allele2)



#----------------------------------------------------------
#----------------------------------------------------------
main ()



