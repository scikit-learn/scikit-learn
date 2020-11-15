#!/bin/bash

set -e
set -x

if [[ "$TRAVIS_CPU_ARCH" != "arm64" ]]; then
    PYTEST="pytest -n $CI_CPU_COUNT" make test-doc
fi
