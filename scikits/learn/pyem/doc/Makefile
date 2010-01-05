# Last Change: Mon May 28 10:00 AM 2007 J

# This makefile is used to build the pdf from the rest file and inlined code
# from python examples

py2tex	= PYTHONPATH=/home/david/local/lib/python2.4/site-packages pygmentize -l python -f tex
rst2tex	= PYTHONPATH=/home/david/local/lib/python2.4/site-packages rst2newlatex.py \
		  --stylesheet-path base.tex --user-stylesheet user.tex 

pytexfiles	= pyem.tex basic_example1.tex basic_example2.tex basic_example3.tex

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

basic_example1.tex: examples/basic_example1.py
	$(py2tex) $< > $@

basic_example2.tex: examples/basic_example2.py
	$(py2tex) $< > $@

basic_example3.tex: examples/basic_example3.py
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
