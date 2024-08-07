#!/bin/bash

set -e
set -x

python -c "import joblib; print(f'Number of cores (physical): \
{joblib.cpu_count()} ({joblib.cpu_count(only_physical_cores=True)})')"

FREE_THREADED_BUILD="$(python -c"import sysconfig; print(bool(sysconfig.get_config_var('Py_GIL_DISABLED')))")"
if [[ $FREE_THREADED_BUILD == "True" ]]; then
    # TODO: delete when importing numpy no longer enables the GIL
    # setting to zero ensures the GIL is disabled while running the
    # tests under free-threaded python
    export PYTHON_GIL=0
fi

# Upgrade numpy to nightly build:
pip install --pre --extra-index https://pypi.anaconda.org/scientific-python-nightly-wheels/simple scikit-learn

# Test that there are no links to system libraries in the
# threadpoolctl output section of the show_versions output:
python -c "import sklearn; sklearn.show_versions()"

# Print information about the machine (Type and generation of CPU, number of
# cores...)
system_profiler SPHardwareDataType

# Dump a core file in case of low level crashes.
echo "Preparing to dump core files"
ulimit -c unlimited
sudo mkdir -p /cores
sudo chmod 1777 /cores

echo "Adding entitlements to python"
/usr/libexec/PlistBuddy -c "Add :com.apple.security.get-task-allow bool true" python.entitlements
/usr/bin/codesign -s - -f --entitlements python.entitlements `which python`

echo "Triggering a segfault manually to check that extracting backtraces works"
python -X faulthandler -c "import ctypes; ctypes.string_at(0)" || (find /cores -name "core.*" -exec lldb -c {} --batch -o "thread backtrace all" -o "quit" \; && exit 0)

echo "Removing the core file from the debug run"
rm -rf /cores/core.*

echo "Running the tests with lldb backtrace reporint on failure"
if pip show -qq pytest-xdist; then
    XDIST_WORKERS=$(python -c "import joblib; print(joblib.cpu_count(only_physical_cores=True))")
    python -X faulthandler -m pytest --pyargs sklearn -n $XDIST_WORKERS || (find /cores -name "core.*" -exec lldb -c {} --batch -o "thread backtrace all" -o "quit" \; && exit 1)
else
    python -X faulthandler -m pytest --pyargs sklearn || (find /cores -name "core.*" -exec lldb -c {} --batch -o "thread backtrace all" -o "quit" \; && exit 1)
fi
