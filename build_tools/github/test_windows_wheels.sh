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
    docker exec -e SKLEARN_SKIP_NETWORK_TESTS=1 \
                -e OMP_NUM_THREADS=2 \
                -e OPENBLAS_NUM_THREADS=2 \
                pytest --pyargs sklearn

    # Test that there are no links to system libraries
    docker exec --rm scikit-learn/minimal-windows \
                python -m threadpoolctl -i sklearn
fi
