#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1
BITNESS=$2

if [[ "$PYTHON_VERSION" == "36" || "$BITNESS" == "32" ]]; then
    # Python 3.6 and 32-bit architectures are not supported
    # by the official Docker images: Tests will just be run
    # on the host (instead of the minimal Docker container).
    exit 0
fi

if [[ "$PYTHON_VERSION" == "37" ]]; then
    WHEEL_NAME=cp"$PYTHON_VERSION"m-win_amd64.whl
else
    # Different name for Python 3.8 and newer
    WHEEL_NAME=cp"$PYTHON_VERSION"-win_amd64.whl
fi

WHEEL_NAME=scikit_learn-$SCIKIT_LEARN_VERSION-cp$PYTHON_VERSION-$WHEEL_NAME

# Find the repaired wheel because there is no a way
# to access to the path in a straightforward manner
WHEEL_PATH="$HOME/AppData/Local/Temp"
WHEEL_PATH=$(find $WHEEL_PATH -type d -name "repaired_wheel")
WHEEL_PATH=$(realpath $WHEEL_PATH)
WHEEL_PATH="$WHEEL_PATH/$WHEEL_NAME"

cp $WHEEL_PATH $WHEEL_NAME

# Dot the Python version for identyfing the base Docker image
PYTHON_VERSION=$(echo ${PYTHON_VERSION:0:1}.${PYTHON_VERSION:1:2})

# Build a minimal Windows Docker image for testing the wheels
docker build --build-arg PYTHON_VERSION=$PYTHON_VERSION \
             --build-arg WHEEL_NAME=$WHEEL_NAME \
             --build-arg CIBW_TEST_REQUIRES="$CIBW_TEST_REQUIRES" \
             -f build_tools/github/Windows \
             -t scikit-learn/minimal-windows .
