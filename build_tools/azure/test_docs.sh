#!/bin/bash

set -ex

source build_tools/shared.sh
activate_environment

scipy_doctest_installed=$(python -c 'import scipy_doctest' && echo "True" || echo "False")
if [[ "$scipy_doctest_installed" == "True" ]]; then
    doc_rst_files=$(find $PWD/doc -name '*.rst' | sort)
    # Changing dir, as we do in build_tools/azure/test_script.sh, avoids an
    # error when importing sklearn. Not sure why this happens ... I am going to
    # wild guess that it has something to do with the bespoke way we set up
    # conda with putting conda in the PATH and source activate, rather than
    # source <conda_root>/etc/profile.d/conda.sh + conda activate.
    cd $TEST_DIR
    # with scipy-doctest, --doctest-modules only runs doctests (in contrary to
    # vanilla pytest where it runs doctests on top of normal tests)
    python -m pytest --doctest-modules --pyargs sklearn
    python -m pytest --doctest-modules $doc_rst_files
fi
