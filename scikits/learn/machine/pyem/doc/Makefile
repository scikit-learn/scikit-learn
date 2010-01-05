# Last Change: Mon Jul 02 09:00 PM 2007 J

# This makefile is used to build the pdf from the rest file and inlined code
# from python examples

py2tex	= PYTHONPATH=/home/david/local/lib/python2.4/site-packages pygmentize -l python -f tex
rst2tex	= PYTHONPATH=/home/david/local/lib/python2.4/site-packages rst2newlatex.py \
		  --stylesheet-path base.tex --user-stylesheet user.tex 

pytexfiles	= pyem.tex basic_example1.tex basic_example2.tex basic_example3.tex pdfestimation.tex discriminant_analysis.tex

SOURCEPATH	= $(PWD)

EXTTOCLEAN=.chk .dvi .log .aux .bbl .blg .blig .ilg .toc .lof .lot .idx .ind .out .bak .ps .pdf .bm

tutorial.pdf: pyem.pdf
	mv $< $@

pyem.pdf: $(pytexfiles)
	pdflatex $<
	pdflatex $<
	pdflatex $<

pyem.tex: index.txt
	$(rst2tex) $< > $@

basic_example1.tex: ../examples/basic_example1.py
	$(py2tex) $< > $@

basic_example2.tex: ../examples/basic_example2.py
	$(py2tex) $< > $@

basic_example3.tex: ../examples/basic_example3.py
	$(py2tex) $< > $@

pdfestimation.tex: ../examples/pdfestimation.py
	$(py2tex) $< > $@

discriminant_analysis.tex: ../examples/discriminant_analysis.py
	$(py2tex) $< > $@

clean:
	for i in $(pytexfiles); do \
		rm -f `echo $$i`; \
	done;
	for i in $(SOURCEPATH); do \
		for j in $(EXTTOCLEAN); do \
			rm -f  `echo $$i/*$$j`; \
		done; \
	done;
