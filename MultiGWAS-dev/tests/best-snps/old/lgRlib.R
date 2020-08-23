
options (width=300, warn=0)
#----------------------------------------------------------
# Util to print head of data
# Fast head, for debug only
#----------------------------------------------------------
hd <- function (data, n=10,m=10) {
	message (deparse (substitute (data)),":")
	if (is.null (dim (data)))
		print (data [1:10])
	else if (ncol (data) < 10) 
		print (data[1:n,])
	else if (nrow (data) < 10)
		print (data[,1:m])
	else 
		print (data [1:n, 1:m])
}
