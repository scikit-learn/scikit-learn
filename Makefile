# simple makefile to simplify repetitive build env management tasks under posix

# caution: testing won't work on windows, see README

PYTHON ?= python
CYTHON ?= cython
PYTEST ?= pytest

# skip doctests on 32bit python
BITS := $(shell python -c 'import struct; print(8 * struct.calcsize("P"))')

all: clean inplace test

clean:
	$(PYTHON) setup.py clean
	rm -rf dist

in: inplace # just a shortcut
inplace:
	$(PYTHON) setup.py build_ext -i

dev-meson:
	pip install --verbose --no-build-isolation --editable . --check-build-dependencies --config-settings editable-verbose=true

clean-meson:
	pip uninstall -y scikit-learn

test-code: in
	$(PYTEST) --showlocals -v sklearn --durations=20
test-sphinxext:
	$(PYTEST) --showlocals -v doc/sphinxext/
test-doc:
ifeq ($(BITS),64)
	$(PYTEST) $(shell find doc -name '*.rst' | sort)
endif
test-code-parallel: in
	$(PYTEST) -n auto --showlocals -v sklearn --durations=20

test-coverage:
	rm -rf coverage .coverage
	$(PYTEST) sklearn --showlocals -v --cov=sklearn --cov-report=html:coverage
test-coverage-parallel:
	rm -rf coverage .coverage .coverage.*
	$(PYTEST) sklearn -n auto --showlocals -v --cov=sklearn --cov-report=html:coverage

test: test-code test-sphinxext test-doc

cython:
	python setup.py build_src

doc: inplace
	$(MAKE) -C doc html

doc-noplot: inplace
	$(MAKE) -C doc html-noplot
