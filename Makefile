# simple makefile to simplify repetitive build env management tasks under posix

PYTHON ?= python

clean:
	$(PYTHON) setup.py clean
	rm -rf dist

inplace:
	$(PYTHON) setup.py build_ext -i

dev-meson:
	pip install --verbose --no-build-isolation --editable . --check-build-dependencies --config-settings editable-verbose=true

clean-meson:
	pip uninstall -y scikit-learn
