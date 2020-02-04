#!/usr/bin/Rscript

library (stringr)
library (dplyr)
options (width=300)
#options(scipen=999)
nROWS = 20

model="Naive"

files =  list.files(".", pattern="^(.*(Naive).*(tbl)[^$]*)$")

files = c("out-Gwasp4-Naive-significativeQTLs.tbl", "out-Gwasp2-Naive-significativeQTLs.tbl", 
		  "out-Plink-Naive-assoc.linear.adjusted.tbl", "out-Tassel-Naive-GLM_Stats_geno+pheno.tbl")
#files = c("out-Gwasp4-Naive-significativeQTLs.tbl", 
summTable = data.frame ()
pThreshold = round (1/(578*2068), 10)
for (f in files) {
	data = read.table (file=f, header=T)
	print (f)
	if (str_detect(f, "Gwasp")) {
		if (str_detect(f, "Gwasp4")) tool = "Gwasp4" else tool = "Gwasp2"
		snps	= data$Marker
		pVal	= round (10^(-data$Score),10)
		chrom   = data$Chrom
		pos	 = data$Position
		signf    = data$Score >= data$Threshold
		pscores = data$Score
		tscores = round (-log10 (pThreshold), 4)
	}else if (str_detect (f, "Plink")) {
		data = data [1:nROWS,]
		tool   = "Plink"
		snps   = data$SNP
		pVal   = data$UNADJ
		chrom  = data$CHR
		pos	= NA
		pscores  = round (-log10 (pVal), 4)
		tscores = round (-log10 (pThreshold), 4)
		signf   = pscores >= tscores
		
	}else if (str_detect (f, "Tassel")) {
		data = data %>% top_n (-1*nROWS, p)
		tool   = "Tassel"
		snps   = data$Marker
		pVal   = data$p
		chrom  = data$Chr
		pos	= data$Pos
		pscores  = round (-log10 (pVal), 4)
		tscores = round (-log10 (pThreshold), 4)
		signf    = pscores >= tscores
	}

	df = data.frame (TOOL=tool, MODEL=model, CHR=chrom, POS=pos, SNPs=snps, P = pVal, pTHR=pThreshold, pSCORE=pscores, tSCORE=tscores, Signf=signf )
	df = df %>% distinct (SNPs, .keep_all=T)
	#df = data.frame (TOOL=tool, MODEL=model, SNPs=snps, P = pVal, CHR=chrom, POS=pos,pTHR=pThreshold)
	summTable = rbind (summTable, df)
	#summTableSorted = summTable %>% arrange (SNPs, P)
	summTableSorted = summTable %>% add_count (SNPs, sort=T, name="N1") %>% arrange (desc(N1))
	summTableSorted = summTableSorted %>% add_count (SNPs, Signf, sort=T, name="Ns") %>% arrange (desc(Ns))
	write.table (file="summary-gwas.tbl", summTable, row.names=F,quote=F, sep="\t")
	write.table (file="summary-gwas-sorted.tbl", summTableSorted, row.names=F,quote=F, sep="\t")
}
