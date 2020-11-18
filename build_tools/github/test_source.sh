#!/bin/bash

set -e
set -x

python -m venv test_env
source test_env/bin/activate

# Run the tests on the installed source distribution
mkdir tmp_for_test
cd tmp_for_test

python -m pip install pytest pandas
python -m pip install dist/*.tar.gz

pytest --pyargs sklearn
