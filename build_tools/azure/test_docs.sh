#!/bin/bash

set -ex

source build_tools/shared.sh
activate_environment
show_installed_libraries

conda info -e
which python
which pytest

echo $TEST_DIR
mkdir -p $TEST_DIR
cp setup.cfg $TEST_DIR
cd $TEST_DIR

python -c 'import sklearn; print(sklearn.__file__)'

scipy_doctest_installed=$(python -c 'import scipy_doctest' && echo "True" || echo "False")
if [[ "$scipy_doctest_installed" == "True" ]]; then
    # With scipy-doctests --doctest-modules only run doctests (in contrary to what happens with vanilla pytest)
    pytest --doctest-modules --pyargs sklearn || echo "failure to run doctests with pyargs"
    pytest --doctest-modules sklearn || echo "failure to run doctests"
    python -m pytest --pyargs sklearn
    pytest $(find doc -name '*.rst' | sort)
fi
