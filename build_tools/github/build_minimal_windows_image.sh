#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1
BITNESS=$2

if [[ "$BITNESS" == "32" ]]; then
    # 32-bit architectures are not supported
    # by the official Docker images: Tests will just be run
    # on the host (instead of the minimal Docker container).
    # the last scipy wheel available for win32 is 1.9.1
    python -m pip install 'scipy<=1.9.1'
    exit 0
fi

TEMP_FOLDER="$HOME/AppData/Local/Temp"
WHEEL_PATH=$(ls -d $TEMP_FOLDER/**/*/repaired_wheel/*)
WHEEL_NAME=$(basename $WHEEL_PATH)

cp $WHEEL_PATH $WHEEL_NAME

# Dot the Python version for identyfing the base Docker image
PYTHON_VERSION=$(echo ${PYTHON_VERSION:0:1}.${PYTHON_VERSION:1:2})

# TODO: Remove when 3.11 images will be available for
# windows (for now the docker image is tagged as 3.11-rc)
# TODO: Remove when 3.11 images will be available for
# scipy (for now only available from nightly anaconda index)
if [[ "$PYTHON_VERSION" == "3.11" ]]; then
    PYTHON_VERSION=$(echo ${PYTHON_VERSION}-rc)
fi

# Build a minimal Windows Docker image for testing the wheels
docker build --build-arg PYTHON_VERSION=$PYTHON_VERSION \
             --build-arg WHEEL_NAME=$WHEEL_NAME \
             --build-arg CONFTEST_NAME=$CONFTEST_NAME \
             --build-arg CIBW_TEST_REQUIRES="$CIBW_TEST_REQUIRES" \
             --build-arg PIP_EXTRA_INDEX_URL="$PIP_EXTRA_INDEX_URL" \
             -f build_tools/github/Windows \
             -t scikit-learn/minimal-windows .
