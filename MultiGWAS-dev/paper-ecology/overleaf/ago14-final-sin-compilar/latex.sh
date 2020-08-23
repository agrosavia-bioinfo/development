DOC=MultiGWASv14
#lyx -e latex $DOC.lyx
pdflatex $DOC.tex
bibtex $DOC
pdflatex $DOC
pdflatex $DOC
#pdf $DOC.pdf
