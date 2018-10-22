#!/bin/bash

set -e

if [[ "$SKIP_TESTS" != "true" ]]; then
    set -x
    make test-doc
fi
