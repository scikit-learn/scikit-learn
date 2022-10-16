#!/bin/bash

# This script is meant to be called by the "install" step
# defined in the ".travis.yml" file. In particular, it is
# important that we call to the right installation script.

if [[ $BUILD_WHEEL == true ]]; then
    source build_tools/travis/install_wheels.sh || travis_terminate 1
else
    source build_tools/travis/install_main.sh || travis_terminate 1
fi
