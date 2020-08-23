#!/usr/bin/Rscript

source ("lglib01.R")
suppressMessages (library ("RColorBrewer"))  # For chord diagrams


createModelsReport <- function (tool) {
	inputDir = "."
	files =  list.files(inputDir, pattern=sprintf ("(^(model).*(%s).*[.](csv))", tool), full.names=T)
	groups = c("N10", "N50", "N100", "N200")

	dfSumm = data.frame ()
	i = 1
	for (f in files) {
		print (f)
		dfUnsorted = read.csv (f)
		df = dfUnsorted [order (dfUnsorted$MODEL,decreasing=T),]

		if (i==1) {
			print (unlist (dfSumm$MODEL))
			dfSumm = data.frame (MODEL=df$MODEL, df$score)
			hd (dfSumm)
		}else {
			models = unlist (dfSumm$MODEL)
			print (models)
			dfSumm = cbind (dfSumm, df [match (models, df$MODEL),"score"])
			hd (dfSumm)
		}
		i = i+1
	}
	outFile = paste0 ("out-best-models-",tool)
	colnames (dfSumm) = c("MODEL", groups)
	write.csv (dfSumm, paste0 (outFile, ".csv"), row.names=F, quote=F)

	pdf (paste0 (outFile, ".pdf"), width=9, height=7)
		# Expand right side of clipping rect to make room for the legend
		par(xpd=T, mar=par()$mar+c(0,0,0,5), cex=1.0)

		
		N = nrow (dfSumm)
		COLORS = brewer.pal (n=N, name="RdBu")
		#COLORS <- colorRampPalette(c("blue", "green", "red"))(n = N)
		mat = as.matrix (dfSumm [,-1])
		barplot (mat, beside=T, col=COLORS, 
				 main=paste0(tool, " models ", " trait total_yield"), ylab="Scores", xlab="Samples")
		legend("right", inset=c(-0.2,0), title="Models", 
			   legend=dfSumm$MODEL, cex=1.0, bty="n", fill=COLORS)
	dev.off()
}

tools = c ("GWASpoly", "TASSEL")

for (t in tools) 
	createModelsReport (t)
