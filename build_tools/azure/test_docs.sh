#!/bin/bash

set -e

if [[ "$DISTRIB" == "conda" ]]; then
    source activate $VIRTUALENV
fi

make test-doc
