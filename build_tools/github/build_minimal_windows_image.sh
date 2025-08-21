#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1

FREE_THREADED_BUILD="$(python -c"import sysconfig; print(bool(sysconfig.get_config_var('Py_GIL_DISABLED')))")"

if [[ $FREE_THREADED_BUILD == "False" ]]; then
    # Prepare a minimal Windows environment without any developer runtime libraries
    # installed to check that the scikit-learn wheel does not implicitly rely on
    # external DLLs when running the tests.
    TEMP_FOLDER="$HOME/AppData/Local/Temp"
    WHEEL_PATH=$(ls -d $TEMP_FOLDER/**/*/repaired_wheel/*)
    WHEEL_NAME=$(basename $WHEEL_PATH)

    cp $WHEEL_PATH $WHEEL_NAME

    # Dot the Python version for identifying the base Docker image
    PYTHON_DOCKER_IMAGE_PART=$(echo ${PYTHON_VERSION:0:1}.${PYTHON_VERSION:1:2})

    if [[ "$CIBW_PRERELEASE_PYTHONS" =~ [tT]rue ]]; then
        PYTHON_DOCKER_IMAGE_PART="${PYTHON_DOCKER_IMAGE_PART}-rc"
    fi

    # Temporary work-around to avoid a loky issue on Windows >= 3.13.7
    if [[ "$PYTHON_DOCKER_IMAGE_PART" == "3.13" ]]; then
        PYTHON_DOCKER_IMAGE_PART="3.13.6"
    fi

    # We could have all of the following logic in a Dockerfile but it's a lot
    # easier to do it in bash rather than figure out how to do it in Powershell
    # inside the Dockerfile ...
    DOCKER_IMAGE="winamd64/python:${PYTHON_DOCKER_IMAGE_PART}-windowsservercore"
    MNT_FOLDER="C:/mnt"
    CONTAINER_ID=$(docker run -it -v "$(cygpath -w $PWD):$MNT_FOLDER" -d $DOCKER_IMAGE)

    function exec_inside_container() {
        docker exec $CONTAINER_ID powershell -Command $1
    }

    exec_inside_container "python -m pip install $MNT_FOLDER/$WHEEL_NAME"
    exec_inside_container "python -m pip install $CIBW_TEST_REQUIRES"

    # Save container state to scikit-learn/minimal-windows image. On Windows the
    # container needs to be stopped first.
    docker stop $CONTAINER_ID
    docker commit $CONTAINER_ID scikit-learn/minimal-windows
else
    # This is too cumbersome to use a Docker image in the free-threaded case
    # TODO When pandas has a release with a Windows free-threaded wheel we can
    # replace the next line with
    # python -m pip install CIBW_TEST_REQUIRES
    python -m pip install pytest
fi
