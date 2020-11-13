#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1
BITNESS=$2

if [[ "$PYTHON_VERSION" == "36" || "$BITNESS" == "32" ]]; then
    # Python 3.6 and 32-bit architectures are not
    # yet supported by the official Docker images
    exit 0
fi

if [[ "$PYTHON_VERSION" == "37" ]]; then
    IDENTIFIER=cp"$PYTHON_VERSION"m-win_amd64.whl
else
    # Different identifier for Python 3.8 and newer
    IDENTIFIER=cp"$PYTHON_VERSION"-win_amd64.whl
fi

WHEEL=scikit_learn-$SCIKIT_LEARN_VERSION-cp$PYTHON_VERSION-$IDENTIFIER
PYTHON_VERSION=$(echo ${PYTHON_VERSION:0:1}.${PYTHON_VERSION:1:2})

docker build --build-arg PYTHON_VERSION=$PYTHON_VERSION \
             --build-arg WHEEL=$WHEEL \
             --build-arg CIBW_TEST_REQUIRES=$CIBW_TEST_REQUIRES
             -t scikit-learn/minimal-windows \
             -f build_tools/github/Windows .
