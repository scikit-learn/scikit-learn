#!/bin/bash

set -e
set -x

FREE_THREADED_BUILD="$(python -c"import sysconfig; print(bool(sysconfig.get_config_var('Py_GIL_DISABLED')))")"
PY_VERSION=$(python -c 'import sys; print(f"{sys.version_info.major}{sys.version_info.minor}")')

if [[ $FREE_THREADED_BUILD == "True" ]]; then
    # TODO: remove when numpy, scipy and pandas have releases with free-threaded wheels
    python -m pip install --pre --extra-index https://pypi.anaconda.org/scientific-python-nightly-wheels/simple numpy scipy pandas --only-binary :all:

elif [[ "$PY_VERSION" == "313" ]]; then
    # TODO: remove when pandas has a release with python 3.13 wheels
    # First install numpy release
    python -m pip install numpy --only-binary :all:
    # Then install pandas-dev
    python -m pip install --pre --extra-index https://pypi.anaconda.org/scientific-python-nightly-wheels/simple pandas --only-binary :all:
fi
