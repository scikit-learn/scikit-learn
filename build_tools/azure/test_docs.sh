#!/bin/bash

set -e

# Defines the show_installed_libraries and activate_environment functions.
source build_tools/shared.sh

activate_environment
make test-doc
