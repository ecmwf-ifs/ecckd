# If you have pdflatex installed, type "make" to create the
# documentation, "make dist-clean" to delete it, "make clean" to just
# remote the latex intermediate files

BASE = ecckd_documentation

documentation: $(BASE).pdf

$(BASE).pdf: $(BASE).tex flowchart.eps
	pdflatex $(BASE).tex
	pdflatex $(BASE).tex
	pdflatex $(BASE).tex

clean:
	rm -f $(BASE).aux $(BASE).log $(BASE).out

dist-clean: clean
	rm -f $(BASE).pdf

.PHONY: documentation clean dist-clean $(BASE).pdf
