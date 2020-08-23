#!/usr/bin/Rscript
library (tableHTML)
library (dplyr)

df = read.table ("out-multiGWAS-inputParameters.tbl", header=T, sep="\t")
nCols = ncol (df)
nRows = nrow (df)

df %>%
  tableHTML(rownames=F,widths=c(360, 200)) %>% 
  add_css_row(css = list('background-color', '#f2f2f2'),
              rows = odd(1:(nRows+1))) %>%
  add_css_row(css = list('background-color', '#e6f0ff'),
              rows = even(1:(nRows+1)))
