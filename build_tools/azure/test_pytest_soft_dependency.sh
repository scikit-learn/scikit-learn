#!/bin/bash

set -e

if [[ "$DISTRIB" == "conda" ]] || [[ "$DISTRIB" == "scipy-dev" ]]; then
    export PATH=$HOME/miniconda3/bin:$PATH
    source activate $VIRTUALENV
elif [[ "$DISTRIB" == "ubuntu" ]]; then
    source $VIRTUALENV/bin/activate
fi

if [[ "$CHECK_PYTEST_SOFT_DEPENDENCY" == "true" ]]; then
    echo "Hello"
fi
