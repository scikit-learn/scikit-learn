#!/bin/bash

set -ex

source build_tools/shared.sh
activate_environment

pwd
which python
which pytest

pytest $(find doc -name '*.rst' | sort)
python -m pytest $(find doc -name '*.rst' | sort)
