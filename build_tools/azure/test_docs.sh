#!/bin/bash

set -ex

source build_tools/shared.sh
activate_environment
show_installed_libraries

conda info -e
which python
which pytest

# Changing dir, as we do in test_script.sh, avoids the fact that we are not
# able to import sklearn. Not sure why this happens ... I am going to wild
# guess that it has something to do with the bespoke way we set up conda with
# putting conda in the PATH and source activate, rather than source
# <conda_root>/etc/profile.d/conda.sh + conda activate.
cd $TEST_DIR

scipy_doctest_installed=$(python -c 'import scipy_doctest' && echo "True" || echo "False")
if [[ "$scipy_doctest_installed" == "True" ]]; then
    python -m pytest --doctest-modules --pyargs sklearn
    python -m pytest $(find doc -name '*.rst' | sort)
fi
