#!/bin/bash

set -e
set -x

FREE_THREADED_BUILD="$(python -c"import sysconfig; print(bool(sysconfig.get_config_var('Py_GIL_DISABLED')))")"
PY_VERSION=$(python -c 'import sys; print(f"{sys.version_info.major}{sys.version_info.minor}")')

# TODO: remove when scipy has a release with free-threaded wheels
if [[ $FREE_THREADED_BUILD == "True" ]]; then
    python -m pip install numpy pandas
    python -m pip install --pre --extra-index https://pypi.anaconda.org/scientific-python-nightly-wheels/simple scipy --only-binary :all:
fi
