# simple makefile to simplify repetitive build env management tasks under posix

PYTHON ?= python

clean-setuptools:
	$(PYTHON) setup.py clean
	rm -rf dist

inplace-setuptools:
	$(PYTHON) setup.py build_ext -i

dev: dev-meson

dev-meson:
	pip install --verbose --no-build-isolation --editable . --check-build-dependencies --config-settings editable-verbose=true

clean-meson:
	pip uninstall -y scikit-learn

clean: clean-meson
