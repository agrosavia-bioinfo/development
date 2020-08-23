
source ("lglib01.R")
library (scales)

mat = read.csv ("x-genoxfeno_3.csv")
hd1 (mat)
ls = unlist (mat)
lsr = rescale (ls, to=c(1,100))
m = matrix (lsr, nrow=nrow(mat), byrow=T)
hd1 (m)
