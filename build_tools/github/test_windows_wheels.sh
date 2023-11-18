#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1

docker container run \
    --rm scikit-learn/minimal-windows \
    powershell -Command "python -c 'import sklearn; sklearn.show_versions()'"

docker container run \
    -e SKLEARN_SKIP_NETWORK_TESTS=1 \
    --rm scikit-learn/minimal-windows \
    powershell -Command "pytest --pyargs sklearn"
