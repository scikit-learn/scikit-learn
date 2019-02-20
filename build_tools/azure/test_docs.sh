#!/bin/bash

set -ex

UNAMESTR=`uname`
if [[ "$UNAMESTR" == "Linux" ]]; then
    source testvenv/bin/activate
fi

# make test-doc
echo "Test Docs"
