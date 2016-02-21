# simple makefile to simplify repetitive build env management tasks under posix

# caution: testing won't work on windows, see README

PYTHON ?= python
CYTHON ?= cython
NOSETESTS ?= nosetests
CTAGS ?= ctags

# skip doctests on 32bit python
BITS := $(shell python -c 'import struct; print(8 * struct.calcsize("P"))')

ifeq ($(BITS),32)
  NOSETESTS:=$(NOSETESTS) -c setup32.cfg
endif


all: clean inplace test

clean-ctags:
	rm -f tags

clean: clean-ctags
	$(PYTHON) setup.py clean
	rm -rf dist

in: inplace # just a shortcut
inplace:
	# to avoid errors in 0.15 upgrade
	rm -f sklearn/utils/sparsefuncs*.so
	rm -f sklearn/utils/random*.so
	$(PYTHON) setup.py build_ext -i

test-code: in
	$(NOSETESTS) -s -v sklearn
test-sphinxext:
	$(NOSETESTS) -s -v doc/sphinxext/
test-doc:
ifeq ($(BITS),64)
	$(NOSETESTS) -s -v doc/*.rst doc/modules/ doc/datasets/ \
	doc/developers doc/tutorial/basic doc/tutorial/statistical_inference \
	doc/tutorial/text_analytics
endif

test-coverage:
	rm -rf coverage .coverage
	$(NOSETESTS) -s -v --with-coverage sklearn

test: test-code test-sphinxext test-doc

trailing-spaces:
	find sklearn -name "*.py" -exec perl -pi -e 's/[ \t]*$$//' {} \;

cython:
	python build_tools/cythonize.py sklearn

ctags:
	# make tags for symbol based navigation in emacs and vim
	# Install with: sudo apt-get install exuberant-ctags
	$(CTAGS) --python-kinds=-i -R sklearn

doc: inplace
	$(MAKE) -C doc html

doc-noplot: inplace
	$(MAKE) -C doc html-noplot

code-analysis:
	flake8 sklearn | grep -v __init__ | grep -v external
	pylint -E -i y sklearn/ -d E1103,E0611,E1101
