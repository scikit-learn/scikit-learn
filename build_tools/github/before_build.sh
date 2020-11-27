#!/bin/bash

set -e
set -x

PYTHON_VERSION=$1

if [[ $PYTHON_VERSION == pp36 ]]; then
    # Install the system build tools
    yum -y -q install openblas-devel
fi
