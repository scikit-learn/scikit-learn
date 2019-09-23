#!/bin/bash

set -e

if [[ "$DISTRIB" =~ ^conda.* ]]; then
    source activate $VIRTUALENV
    conda install scikit-image -y
elif [[ "$DISTRIB" == "ubuntu" ]]; then
    source $VIRTUALENV/bin/activate
    pip install scikit-image
fi

make test-doc
