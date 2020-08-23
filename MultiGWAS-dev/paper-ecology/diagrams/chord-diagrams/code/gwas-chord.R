#!/opt/bin/Rscript

library(circlize)
options (width=300)

args = commandArgs (trailingOnly=T)
args = c ("out-multiGWAS-scoresTable-best.scores")

scoresFile  = args [1]
scores      = read.table (file=scoresFile, sep="\t", header=T)
scores [,3] = paste0 ("Chrom", scores[,3])

tbl = scores [,c(1,3,5)]

# Group by TOOL and select 3 SNPs for each one
#tblr  = Reduce (rbind, by (tbl, tbl["TOOL"], head, n=5))
tblr  = tbl
tblm = tblr [,c(2,3)]
chrs = sort (tblm [!duplicated (tblm [,1]), 1])
snps = sort (as.character (tblm [!duplicated (tblm [,2]), 2]))

# Create matrix Chroms X SNPs
n = length (chrs)
m = length (snps)
mat  = as.data.frame (matrix (rep (0,n*m),nrow=n, ncol=m), stringAsFactor=F )
rownames (mat) = chrs
colnames (mat) = snps

# Fill the matrix

dmat = as.data.frame (mat)
for (i in 1:nrow (tblm)) {
	chr = as.character (tblm [i, 1])
	snp = as.character (tblm [i, 2])
	dmat [chr,snp] = dmat [chr,snp] + 1
}

outFile = paste0(strsplit (scoresFile, split="[.]")[[1]][1], "-chordDiagram.pdf")
pdf (file=outFile, width=7, height=7)

grid.col <- setNames(rainbow(length(unlist(dimnames(dmat))), start=.3, end=.1,
							 ), union(rownames(dmat), colnames(dmat)))

dmat = as.matrix (dmat) 
chordDiagram(dmat, annotationTrack = "grid", directional = -1, direction.type = c("diffHeight", "arrows"),
			 link.arr.type = "big.arrow", reduce=0.01,
    preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(dmat))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

#chordDiagram (as.matrix (dmat), 
#              directional = -1, direction.type = "arrows", link.arr.col = "black",  link.arr.length = 0.5, transparency = 0.2, link.sort = TRUE, link.decreasing = TRUE)
dev.off()

write.table (file="matChrSnp.tbl", dmat, quote=F, sep="\t", row.names=T)

