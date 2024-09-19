#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1

TEMP_FOLDER="$HOME/AppData/Local/Temp"
WHEEL_PATH=$(ls -d $TEMP_FOLDER/**/*/repaired_wheel/*)
WHEEL_NAME=$(basename $WHEEL_PATH)

cp $WHEEL_PATH $WHEEL_NAME

# Dot the Python version for identifying the base Docker image
PYTHON_VERSION=$(echo ${PYTHON_VERSION:0:1}.${PYTHON_VERSION:1:2})

if [[ "$CIBW_PRERELEASE_PYTHONS" =~ [tT]rue ]]; then
    PYTHON_VERSION="$PYTHON_VERSION-rc"
fi

# We could have all of the following logic in a Dockerfile but it's a lot
# easier to do it in bash rather than figure out how to do it in Powershell
# inside the Dockerfile ...
DOCKER_IMAGE="winamd64/python:$PYTHON_VERSION-windowsservercore"
MNT_FOLDER="C:\mnt"
container_id=$(docker run -it -v "$(cygpath -w $PWD):$MNT_FOLDER" -d $DOCKER_IMAGE)

function exec_inside_container() {
    docker exec $container_id powershell -Command $1
}

exec_inside_container "python -m pip install $MNT_FOLDER/$WHEEL_NAME"

if [[ "$PYTHON_VERSION" == "3.13" ]]; then
    # TODO: remove when pandas has a release with python 3.13 wheels
    # First install numpy release
    exec_inside_container "python -m pip install numpy"
    # Then install pandas-dev
    exec_inside_container "python -m pip install --pre --extra-index https://pypi.anaconda.org/scientific-python-nightly-wheels/simple pandas --only-binary :all:"
fi

exec_inside_container "python -m pip install $CIBW_TEST_REQUIRES"

# Save container state to scikit-learn/minimal-windows image
docker commit $container_id scikit-learn/minimal-windows
