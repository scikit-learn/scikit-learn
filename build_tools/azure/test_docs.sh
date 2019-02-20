#!/bin/bash

set -ex

if [[ "$DISTRIB" == "ubuntu" ]]; then
    source testvenv/bin/activate
fi

# make test-doc
echo "Test Docs"
