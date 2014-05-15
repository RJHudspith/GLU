## makes the documentation
doxygen:
	doxygen Doxyfile

.PHONY: 
	doxygen

documentation:
	pdflatex Setup_info.tex && pdflatex Setup_info_wflow.tex && pdflatex Gluon_info.tex && pdflatex Smearing_info.tex && pdflatex SUNxU1_info.tex && pdflatex Beginnings.tex && mv -v *.pdf ${DOCPREFIX}/docs/ && rm *.out *.log *.aux *~

.PHONY:
	documentation