#!/usr/bin/Rscript

# Extract and sort the names with integer suffix from a table

#args = commandArgs(trailingOnly = TRUE)
args = "ClusterCall_prediction_CCC-fixed.csv"
print (args)
filename    = args [1]
outFilename = sprintf ("%s-names", strsplit (filename, "[.]")[[1]][1]) 
print (sprintf ("\nExtracting names from %s to %s", filename, outFilename))

#filename = "berdugo-genotype-s003.csv"

# Load and get names
genoB.all = read.table (filename, header=T, sep=",")
dim (genoB.all)
nms = names (genoB.all)


# Filter to specific columns
nms.andall = nms [grep ("And", nms)]
#nms.and = nms.andall [-grep ("[.]", nms.andall)]
nms.and = nms.andall

# Create ids from suffixes
la = unlist (strsplit(nms.and, "_"))
nms.ids = sprintf ("%03d", strtoi (la [la != "And"]))

# Create matrix ids + names, sort, and write
m =  matrix(c(nms.ids,nms.and), ncol=2)
markers.sorted = m[order (m[,1]),]
write.table (file=sprintf ("%s-ids.tbl", outFilename), markers.sorted,row.names=F, quote=F,col.names=F)
write.table (file=sprintf ("%s.tbl", outFilename), markers.sorted[,2],row.names=F, quote=F,col.names=F)
