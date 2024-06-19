#!/bin/bash

set -e

source build_tools/shared.sh
activate_environment

make test-doc
