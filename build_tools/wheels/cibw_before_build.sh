#!/bin/bash

set -euxo pipefail

PROJECT_DIR="$1"
LICENSE_FILE="$PROJECT_DIR/COPYING"

echo "" >>"$LICENSE_FILE"
echo "----" >>"$LICENSE_FILE"
echo "" >>"$LICENSE_FILE"

if [[ $RUNNER_OS == "Linux" ]]; then
    cat $PROJECT_DIR/build_tools/wheels/LICENSE_linux.txt >>"$LICENSE_FILE"
elif [[ $RUNNER_OS == "macOS" ]]; then
    cat $PROJECT_DIR/build_tools/wheels/LICENSE_macos.txt >>"$LICENSE_FILE"
elif [[ $RUNNER_OS == "Windows" ]]; then
    cat $PROJECT_DIR/build_tools/wheels/LICENSE_windows.txt >>"$LICENSE_FILE"
fi

PYTHON_VERSION=$(python -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
IS_PYTHON_FREE_THREADED=$(python -c 'import sysconfig; print(sysconfig.get_config_var("Py_GIL_DISABLED"))')

# Use abi3 wheels for non vanilla Python >= 3.12
if [[ $PYTHON_VERSION != "3.11" && IS_PYTHON_FREE_THREADED == "0" ]]; then
    echo "CIBW_CONFIG_SETTINGS='setup-args=-Dpython.allow_limited_api=true" >> "$GITHUB_ENV"
fi
