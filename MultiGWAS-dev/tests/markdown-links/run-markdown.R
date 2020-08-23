#!/usr/bin/Rscript
.libPaths ("/home/lg/agrosavia/GWAS/MultiGWAS-dev/opt/Rlibs")

rmarkdown::render ("gwas-markdown.Rmd", output_file="out.html", output_format="html_document", 
				   #output_options=list(self_contained=F),
				   params=list (workingDir=getwd(), reportTitle="Naive report", nBest=5))

