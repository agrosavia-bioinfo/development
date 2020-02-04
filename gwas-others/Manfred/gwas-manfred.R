#!/usr/bin/Rscript

suppressMessages (library (dplyr))
options (width=300, stringsAsFactors=T)
args = commandArgs(trailingOnly=T)

args = c ("manfred-mca2012.tbl", "samples-filtered-names.tbl")
traitsFile  = args [1]
samplesFile = args [2]


data2012 = read.csv (file=traitsFile, row.names=1)
samples  = as.character (read.table (samplesFile,header=T)[,1])

data = data2012 [row.names (data2012) %in% samples, -c(1,2,3)]
print (data [1:20,1:30])

#data = data %>% select (

