#Command to run SMARTPCA:
#To use this script, make sure that program SMARTPCA  is executable and saved with the name "smartpca" in the same folder from where you run fcGENE.
#If this is not the case, you should give its  path and name while excuting above commands.
#If necessary, you can also edit (add/delete) the necessary command options in file: fcgene_out_smartpca.par
system("./smartpca -p fcgene_out_smartpca.par>fcgene_out_smartpca_screenOutput.txt")
#R function to make pca plot: 
pca_plot<-function(
		evecFileName, 	#name of *.evec file obtained from eigenstrat
		outFileName	="", #name of plot file
		used_evecs=1:2, #which eigenvectors to plot.  Use 1:2 for top two eigenvectors. 
		file_type="pdf" #type of  file to be plotted.You can also write  "png".
	 	){
 		if(outFileName==""){
 			outFileName=substring(evecFileName,1,(nchar(evecFileName)-nchar(".evec")))
	 		plotName<-paste(outFileName, ".",file_type,sep="")
		}
	 	evec<-read.table(evecFileName)
	 	ncols<-dim(evec)[2]
		evec<-evec[c(used_evecs+1,ncols)]
	 	colnames(evec)<-c("pca1","pca2","groupName")
	 	groupNames<-evec$groupName
	 	color<-1:length(unique(groupNames))
	 	evec$color<-color[groupNames]
	 	evec$pch<-substring(groupNames,1,5)
		if(file_type=="pdf")
	 		pdf(plotName)
		else if(file_type=="png")
	 	 	 png(plotName)
		xvalue<-as.numeric(evec$pca1)
	 	yvalue<- as.numeric(evec$pca2)
	 	plot(xvalue,yvalue , pch="",xlab="first_component",ylab="second_component")
	 	text(xvalue,yvalue, labels=evec$pch,col=evec$color)
		dev.off()
	 	cat("pca plot can be found in file:\t",plotName,"\n")
}

evecFileName<-"fcgene_out.evec"
pca_plot(evecFileName)
#To write file"fcgene_out.pca" necessary for eigenstrat: 
evec<-as.matrix(read.table(evecFileName,comment.char="",fill=T,header=F))
evec<-evec[,c(-1,-dim(evec)[2])]
outFileName<-"fcgene_out.pca"
sink(outFileName) #this command is to open file for writing!
cat(dim(evec)[2],"
")#writes number of eigenvectors in the file 
writeLines(evec[1,])# write eigenvalues given in *.evec file 
write.table(evec[-1,],file=outFileName,append=T,quot=F,col.names=F,row.names=F,sep=" " )
sink() #file closing part 
unlist(outFileName)#file closing part
