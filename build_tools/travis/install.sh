#!/bin/bash

# This script is meant to be called by the "install" step
# defined in the ".travis.yml" file. In particular, it is
# important that we call to the right installation script.

set -e

if [[ $BUILD_WHEEL == true ]]; then
    source build_tools/travis/install_wheels.sh
else
    source build_tools/travis/install_main.sh
fi
