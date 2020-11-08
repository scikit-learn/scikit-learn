#!/bin/bash

set -e
set -x

python -m pip install numpy scipy cython
python -m pip install twine
python -m pip install pytest pandas

python setup.py sdist
python -m pip install dist/*.tar.gz
python setup.py build_ext -i

pytest --pyargs sklearn

# Check whether the source distribution will render correctly
twine check dist/*.tar.gz
