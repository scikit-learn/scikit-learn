#!/bin/bash

set -e

source build_tools/shared.sh
activate_environment

python -m pytest $(find doc -name '*.rst' | sort)
