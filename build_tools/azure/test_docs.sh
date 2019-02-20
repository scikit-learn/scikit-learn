#!/bin/bash

set -e

UNAMESTR=`uname`
if [[ "$UNAMESTR" == "Linux" ]]; then
    source $VIRTUALENV_DIR/bin/activate
fi

# make test-doc
echo "Test Docs"
