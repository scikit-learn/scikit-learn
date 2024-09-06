#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1

if [[ "$PYTHON_VERSION" == "313" ]]; then
    # TODO: remove when pandas has a release with python 3.13 wheels
    # First install numpy release
    docker container run \
        --rm scikit-learn/minimal-windows \
        powershell -Command "python -m pip install numpy"
    # Then install pandas-dev
    docker container run \
        --rm scikit-learn/minimal-windows \
        powershell -Command "python -m pip install --pre --extra-index https://pypi.anaconda.org/scientific-python-nightly-wheels/simple pandas --only-binary :all:"
fi

docker container run \
    --rm scikit-learn/minimal-windows \
    powershell -Command "python -m pip install $CIBW_TEST_REQUIRES"

docker container run \
    --rm scikit-learn/minimal-windows \
    powershell -Command "python -m pip install  -c 'import sklearn; sklearn.show_versions()'"

docker container run \
    --rm scikit-learn/minimal-windows \
    powershell -Command "python -c 'import sklearn; sklearn.show_versions()'"

docker container run \
    -e SKLEARN_SKIP_NETWORK_TESTS=1 \
    --rm scikit-learn/minimal-windows \
    powershell -Command "pytest --pyargs sklearn"
