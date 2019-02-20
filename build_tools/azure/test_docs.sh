#!/bin/bash

set -e

if [[ "$DISTRIB" == "conda" ]]; then
    export PATH=$MINICONDA_PATH/bin:$PATH
    source activate $VIRTUALENV
elif [[ "$DISTRIB" == "ubuntu" ]]; then
    source $VIRTUALENV/bin/activate
fi

make test-doc
