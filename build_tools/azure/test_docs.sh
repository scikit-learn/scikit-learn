#!/bin/bash

set -ex

source build_tools/shared.sh
activate_environment

scipy_doctest_installed=$(python -c 'import scipy_doctest' && echo "True" || echo "False")
if [[ "$scipy_doctest_installed" == "True" ]]; then
    # With scipy-doctests --doctest-modules only run doctests (in contrary to what happens with vanilla pytest)
    python -m pytest --doctest-modules --pyargs sklearn
    python -m pytest $(find doc -name '*.rst' | sort)
fi
