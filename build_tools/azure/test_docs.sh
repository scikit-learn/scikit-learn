#!/bin/bash

set -ex

source build_tools/shared.sh
activate_environment

# XXX: for some unknown reason python -m pytest fails here in the CI, can't
# reproduce locally and not worth spending time on this
pytest $(find doc -name '*.rst' | sort)
