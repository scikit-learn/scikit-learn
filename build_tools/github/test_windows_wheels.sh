#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1
BITNESS=$2

if [[ "$PYTHON_VERSION" == "36" || "$BITNESS" == "32" ]]; then
    # For Python 3.6 and 32-bit architecture use the regular
    # test command (outside of the minimal Docker container)
    pytest --pyargs sklearn

    # Test that there are no links to system libraries
    python -m threadpoolctl -i sklearn
else
    docker run --rm scikit-learn/minimal_windows
fi
