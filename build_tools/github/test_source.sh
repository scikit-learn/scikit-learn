#!/bin/bash

set -e
set -x

cd ../../

python -m venv test_env
source test_env/bin/activate

python -m pip install scikit-learn/scikit-learn/dist/*.tar.gz
python -m pip install pytest pandas

# Run the tests on the installed source distribution
mkdir tmp_for_test
cp scikit-learn/scikit-learn/conftest.py tmp_for_test
cd tmp_for_test

pytest --pyargs sklearn
