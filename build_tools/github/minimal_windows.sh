#!/bin/bash

set -e
set -x

# Remove the dot from the Python version for naming the wheel
TRIM_PYTHON_VERSION=$(echo $PYTHON_VERSION | tr -d ".")

if [[ $PYTHON_VERSION == 3.7 ]]; then
    IDENTIFIER=cp"$TRIM_PYTHON_VERSION"m-win_amd64.whl
else
    # Different identifier for Python 3.8 and newer
    IDENTIFIER=cp"$TRIM_PYTHON_VERSION"-win_amd64.whl

WHEEL=scikit_learn-$SCIKIT_LEARN_VERSION-cp$TRIM_PYTHON_VERSION-$IDENTIFIER

docker build --build-arg PYTHON_VERSION=$PYTHON_VERSION \
             --build-arg WHEEL=$WHEEL \
             -t scikit-learn/minimal-windows \
             -f build_tools/github/Windows .

docker run --rm scikit-learn/minimal-windows
