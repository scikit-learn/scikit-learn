#!/bin/bash

set -e
set -x

if [[ "$BUILD_WITH_ICC" == "true" ]]; then
    # the tools in the oneAPI toolkits are configured via environment variables
    # which are also required at runtime.
    source /opt/intel/inteloneapi/setvars.sh
fi

make test-doc
