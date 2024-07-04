# simple makefile to simplify repetitive build env management tasks under posix

PYTHON ?= python
DEFAULT_MESON_BUILD_DIR = build/cp$(shell python -c 'import sys; print(f"{sys.version_info.major}{sys.version_info.minor}")' )

all:
	@echo "Please use 'make <target>' where <target> is one of"
	@echo "  dev                  build scikit-learn with Meson"
	@echo "  clean                clean scikit-learn Meson build. Very rarely needed,"
	@echo "                       one use case is when switching back to setuptools"
	@echo "  dev-setuptools       build scikit-learn with setuptools (deprecated)"
	@echo "  clean-setuptools     clean scikit-learn setuptools build (deprecated)"

.PHONY: all

dev: dev-meson

dev-meson:
	pip install --verbose --no-build-isolation --editable . --check-build-dependencies --config-settings editable-verbose=true

clean: clean-meson

clean-meson:
	pip uninstall -y scikit-learn
	# It seems in some cases (e.g. weird compilation errors when switching
	# between numpy<2 and numpy>=2) removing the folder avoids weird compilation
	# errors. For some reason ninja clean -C $(DEFAULT_MESON_BUILD_DIR) is not
	# enough
	rm -rf $(DEFAULT_MESON_BUILD_DIR)

dev-setuptools:
	$(PYTHON) setup.py build_ext -i

clean-setuptools:
	$(PYTHON) setup.py clean
	rm -rf dist
