SPHINXOPTS=-v

all: clean html show

clean:
	rm -rf _build/*
	rm -rf auto_examples/
	rm -rf auto_examples_README_header/
	rm -rf auto_examples_rst_index/
	rm -rf auto_examples_with_rst/
	rm -rf gen_modules/

html:
	sphinx-build -b html -d _build/doctrees $(SPHINXOPTS) . _build/html

latexpdf:
	sphinx-build -b latex -d _build/doctrees . _build/latex
	@echo "Running LaTeX files through pdflatex..."
	$(MAKE) -C _build/latex all-pdf
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

show:
	@python -c "import webbrowser; webbrowser.open_new_tab('file://$(abspath .)/_build/html/index.html')"
