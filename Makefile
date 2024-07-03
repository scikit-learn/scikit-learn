# simple makefile to simplify repetitive build env management tasks under posix

PYTHON ?= python

all:
	@echo "Please use 'make <target>' where <target> is one of"
	@echo "  dev                  build scikit-learn with Meson"
	@echo "  clean                clean scikit-learn Meson build. Very rarely needed,"
	@echo "                       since meson-python recompiles on import"

.PHONY: all

dev: dev-meson

dev-meson:
	pip install --verbose --no-build-isolation --editable . --check-build-dependencies --config-settings editable-verbose=true

clean: clean-meson

clean-meson:
	pip uninstall -y scikit-learn
