#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1

WHEEL_PATH=$(ls wheelhouse/*.whl)
WHEEL_NAME=$(basename "$WHEEL_PATH")

# Dot the Python version for identifying the base Docker image.
PYTHON_DOCKER_IMAGE_PART=$(echo "${PYTHON_VERSION:0:1}.${PYTHON_VERSION:1:2}")

DOCKER_IMAGE="winamd64/python:${PYTHON_DOCKER_IMAGE_PART}-windowsservercore"
MNT_FOLDER="C:/mnt"
CONTAINER_ID=$(docker run -it -v "$(cygpath -w "$PWD"):$MNT_FOLDER" -d "$DOCKER_IMAGE")

function exec_inside_container() {
    docker exec "$CONTAINER_ID" powershell -NoProfile -Command "$1"
}

exec_inside_container "python -m venv C:/venv"
exec_inside_container "C:/venv/Scripts/python -m pip install $MNT_FOLDER/wheelhouse/$WHEEL_NAME"
exec_inside_container "C:/venv/Scripts/python -c 'import sklearn; sklearn.show_versions()'"
exec_inside_container "C:/venv/Scripts/python -c 'from sklearn.utils._openmp_helpers import _openmp_parallelism_enabled; assert _openmp_parallelism_enabled()'"
