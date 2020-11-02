#!/bin/bash

set -e
set -x

if [[ "$BUILD_WITH_ICC" == "true" ]]; then
    # the tools in the oneAPI toolkits are configured via environment variables
    # which are also required at runtime.
    source /opt/intel/inteloneapi/setvars.sh
fi

if [[ "$TRAVIS_CPU_ARCH" != "arm64" ]]; then
    PYTEST="pytest -n $CI_CPU_COUNT" make test-doc
fi
