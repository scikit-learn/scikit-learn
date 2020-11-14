#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1
BITNESS=$2

if [[ "$PYTHON_VERSION" == "36" || "$BITNESS" == "32" ]]; then
    # For Python 3.6 and 32-bit architecture use the regular
    # test command (outside of the minimal Docker container)
    pytest --pyargs sklearn
    python -m threadpoolctl -i sklearn
else
    docker run -e SKLEARN_SKIP_NETWORK_TESTS=1 \
               --rm scikit-learn/minimal-windows
fi
