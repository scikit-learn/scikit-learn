#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1

# Install the system build tools for PyPy implementation
if [[ $PYTHON_VERSION == pp36 ]]; then
    apt-get -yq install libatlas-base-dev liblapack-dev gfortran libopenblas-dev
fi
