#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1
BITNESS=$2

if [[ "$PYTHON_VERSION" == "36" || "$BITNESS" == "32" ]]; then
    # For Python 3.6 and 32-bit architecture use the regular
    # test command (outside of the minimal Docker container)
    cp $CONFTEST_PATH $CONFTEST_NAME
    pytest --pyargs sklearn
    python -m threadpoolctl -i sklearn
else
    docker container run --name minimal_windows \
                         -d -ti --rm scikit-learn/minimal-windows powershell

    docker exec -e SKLEARN_SKIP_NETWORK_TESTS=1 \
                -e OMP_NUM_THREADS=2 \
                -e OPENBLAS_NUM_THREADS=2 \
                minimal_windows pytest --pyargs sklearn

    docker exec minimal_windows python -m threadpoolctl -i sklearn
fi
