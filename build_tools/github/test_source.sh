#!/bin/bash

set -e
set -x

# Run the tests on the installed source distribution
mkdir test_folder
cp conftest.py test_folder
pushd test_folder
pytest --pyargs sklearn
popd

# Check whether the source distribution will render correctly
twine check dist/*.tar.gz
