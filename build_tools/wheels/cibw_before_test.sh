#!/bin/bash

set -e
set -x

FREE_THREADED_BUILD="$(python -c"import sysconfig; print(bool(sysconfig.get_config_var('Py_GIL_DISABLED')))")"
if [[ $FREE_THREADED_BUILD == "True" ]]; then
    # TODO: remove when numpy, scipy and pandas have releases with free-threaded wheels
    python -m pip install --pre --extra-index https://pypi.anaconda.org/scientific-python-nightly-wheels/simple numpy scipy pandas --only-binary :all:
fi
