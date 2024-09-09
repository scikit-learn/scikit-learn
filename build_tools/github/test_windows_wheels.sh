#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1

container_id=$(docker run -it -d scikit-learn/minimal-windows)

if [[ "$PYTHON_VERSION" == "313" ]]; then
    # TODO: remove when pandas has a release with python 3.13 wheels
    # First install numpy release
    docker exec $container_id \
        powershell -Command "python -m pip install numpy"
    # Then install pandas-dev
    docker exec $container_id \
        powershell -Command "python -m pip install --pre --extra-index https://pypi.anaconda.org/scientific-python-nightly-wheels/simple pandas --only-binary :all:"
fi

docker exec $container_id \
    powershell -Command "python -m pip install $CIBW_TEST_REQUIRES"

docker exec $container_id \
    powershell -Command "python -c 'import sklearn; sklearn.show_versions()'"

docker exec $container_id \
    powershell -Command "python -m pip list" || echo pip list failed

docker exec $container_id \
    -e SKLEARN_SKIP_NETWORK_TESTS=1 \
    powershell -Command "python -m pytest --pyargs sklearn"
