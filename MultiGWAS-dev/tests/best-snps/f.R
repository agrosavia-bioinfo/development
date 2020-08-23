f <- function(x){
		ploidy=4
		x = unlist (x)
		ref   = x[1]
		genos = x[-1]
        y   <- gregexpr(pattern=ref,text=genos,fixed=T)
        ans <- as.integer(lapply(y,function(z){ifelse(z[1]<0,ploidy,ploidy-length(z))}))
        return(ans)
}
