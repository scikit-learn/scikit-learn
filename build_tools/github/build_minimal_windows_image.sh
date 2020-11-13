#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1
BITNESS=$2

if [[ "$PYTHON_VERSION" == "36" || "$BITNESS" == "32" ]]; then
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

IDENTIFIER=scikit_learn-$SCIKIT_LEARN_VERSION-cp$PYTHON_VERSION-$IDENTIFIER

# Find the repaired wheel because there is no a way
# to access to the path in a straightforward manner
BASE_REPAIRED_WHEEL_PATH="C:/Users/RUNNER~1/AppData/Local/Temp"
REPAIRED_WHEEL_PATH=$(find $BASE_PATH -type d -name "repaired_wheel")
REPAIRED_WHEEL_PATH =$(realpath $REPAIRED_WHEEL_PATH)

WHEEL="$PATH/$IDENTIFIER"

# Dot the Python version for identyfing the base Docker image
PYTHON_VERSION=$(echo ${PYTHON_VERSION:0:1}.${PYTHON_VERSION:1:2})

# Build a minimal Windows Docker image for testing the wheels
docker build --build-arg PYTHON_VERSION=$PYTHON_VERSION \
             --build-arg WHEEL="$WHEEL" \
             --build-arg CIBW_TEST_REQUIRES="$CIBW_TEST_REQUIRES" \
             -f build_tools/github/Windows \
             -t scikit-learn/minimal-windows .
