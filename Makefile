# simple makefile to simplify repetitive build env management tasks under posix

PYTHON ?= python

dev: dev-meson

dev-meson:
	pip install --verbose --no-build-isolation --editable . --check-build-dependencies --config-settings editable-verbose=true

clean: clean-meson

clean-meson:
	pip uninstall -y scikit-learn

dev-setuptools:
	$(PYTHON) setup.py build_ext -i

clean-setuptools:
	$(PYTHON) setup.py clean
	rm -rf dist
