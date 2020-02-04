#!/usr/bin/Rscript

options(java.parameters = c("-Xmx1g", "-Xms1g"))
library(rTASSEL)
rTASSEL::startLogger(fullPath = NULL, fileName = NULL)

genoPathHMP <- "mdp_genotype.hmp.txt"
# Load in hapmap file
tasGenoHMP <- rTASSEL::readGenotypeTableFromPath(path=genoPathHMP)
print (tasGenoHMP)

print (">>>>>")
# Load into rTASSEL "GenotypePhenotype" object
tasPheno <- rTASSEL::readPhenotypeFromPath(path="mdp_traits.txt")
print (tasPheno)
