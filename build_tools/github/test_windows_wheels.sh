#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1
BITNESS=$2

if [[ "$PYTHON_VERSION" == "36" || "$BITNESS" == "32" ]]; then
    # Python 3.6 and 32-bit architectures are not
    # yet supported by the official Docker images
    $CIBW_TEST_COMMAND
else
    docker run --rm scikit-learn/minimal_windows
fi
