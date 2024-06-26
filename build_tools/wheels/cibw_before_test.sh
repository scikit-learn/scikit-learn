#!/bin/bash

set -e
set -x

FREE_THREADED_BUILD="$(python -c"import sysconfig; print(bool(sysconfig.get_config_var('Py_GIL_DISABLED')))")"
if [[ $FREE_THREADED_BUILD == "True" ]]; then
    # TODO: remove when numpy and scipy have releases with free-threaded wheels
    python -m pip install --pre --extra-index https://pypi.anaconda.org/scientific-python-nightly-wheels/simple numpy scipy
else
    # There is no pandas free-threaded wheel at the time of writing, so we only
    # install pandas in other builds
    # TODO: adapt when there is a pandas free-threaded wheel
    python -m pip install pandas
fi
