#!/bin/bash

set -e
set -x

python -m pip install numpy scipy cython
python -m pip install twine
python -m pip install pytest pandas

python setup.py sdist
python -m pip install dist/*.tar.gz
