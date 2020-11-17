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
    docker run --name minimal_windows \
               -d -ti --rm scikit-learn/minimal-windows powershell

    # The "-e" option is not working for docker
    # run nor docker exec on Windows containers
    docker exec minimal_windows \$env:SKLEARN_SKIP_NETWORK_TESTS="1"
    docker exec minimal_windows \$env:OMP_NUM_THREADS="2"
    docker exec minimal_windows \$env:OPENBLAS_NUM_THREADS="2"

    docker exec minimal_windows pytest --pyargs sklearn

    # Test that there are no links to system libraries
    docker exec minimal_windows python -m threadpoolctl -i sklearn
fi
