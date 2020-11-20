#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1

if [[ $PYTHON_VERSION == pp36 ]]; then
    # Install the system build tools (ATLAS and LAPLACK) for PyPy implementation
    apt-get -yq install libatlas-base-dev liblapack-dev gfortran libopenblas-dev
fi
