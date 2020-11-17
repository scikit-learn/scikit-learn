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
    docker container run --env SKLEARN_SKIP_NETWORK_TESTS=1 \
                         --env OMP_NUM_THREADS=2 \
                         --env OPENBLAS_NUM_THREADS=2 \
                         --name minimal_windows \
                         --rm scikit-learn/minimal-windows powershell 'Write-Output "$env:SKLEARN_SKIP_NETWORK_TESTS $env:OMP_NUM_THREADS $env:OPENBLAS_NUM_THREADS"'

    # docker exec minimal_windows pytest --pyargs sklearn

    # Test that there are no links to system libraries
    # docker exec minimal_windows python -m threadpoolctl -i sklearn
fi
