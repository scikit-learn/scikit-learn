#!/bin/bash

set -e

if [[ "$DISTRIB" =~ ^conda.* ]]; then
    source activate $VIRTUALENV
    conda install scikit-image -y
elif [[ "$DISTRIB" == "ubuntu" ]]; then
    source $VIRTUALENV/bin/activate
    python -m pip install scikit-image==0.12.3
fi

make test-doc
