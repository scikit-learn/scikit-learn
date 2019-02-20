#!/bin/bash

set -e

UNAMESTR=`uname`
if [[ "$UNAMESTR" == "Linux" ]]; then
    source $VIRTUALENV/bin/activate
fi

make test-doc
