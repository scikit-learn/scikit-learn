#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1
BITNESS=$2

if [[ "$BITNESS" == "32" ]]; then
    # 32-bit architectures use the regular
    # test command (outside of the minimal Docker container)
    cp $CONFTEST_PATH $CONFTEST_NAME
    python -c "import sklearn; sklearn.show_versions()"
    pytest --pyargs sklearn
else
    docker container run \
        --rm scikit-learn/minimal-windows \
        powershell -Command "python -c 'import sklearn; sklearn.show_versions()'"

    docker container run \
        -e SKLEARN_SKIP_NETWORK_TESTS=1 \
        -e OMP_NUM_THREADS=2 \
        -e OPENBLAS_NUM_THREADS=2 \
        --rm scikit-learn/minimal-windows \
        powershell -Command "pytest --pyargs sklearn"
fi
