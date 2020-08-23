#!/usr/bin/Rscript

library(ggplot2)
library(gridExtra)

df = read.table ("out-multiGWAS-inputParameters.tbl", header=T, sep="\t")


mytheme <- ttheme_default(core = list(fg_params = list(hjust=0, x=0.1, fontsize=8)),
                          colhead = list(fg_params = list(fontsize=9, fontface="bold")))

#png("out-test.png")
#tbl<-tableGrob(df, theme=mytheme)
#grid.arrange(plots, table, heights=c(5,1))
#grid.arrange(p)
#dev.off()

library(gridExtra)
png("test.png")
grid.table(df, theme=mytheme)
dev.off()
