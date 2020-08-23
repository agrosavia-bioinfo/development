DOC=MultiGWASv12
#lyx -e latex $DOC.lyx
pdflatex $DOC.tex
bibtex MultiGWASv12-refs
pdflatex $DOC
pdflatex $DOC
pdf $DOC.pdf
