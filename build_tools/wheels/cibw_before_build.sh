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

# Do not use abi3 wheels for Python 3.11 or free-threaded Python
if [[ $PYTHON_VERSION == "3.11" || IS_PYTHON_FREE_THREADED == "1" ]]; then
    echo "CIBW_CONFIG_SETTINGS='setup-args=-Dpython.allow_limited_api=false" >> "$GITHUB_ENV"
else
    # this is necessary for the wheel to be named correctly, is there a better way???
    cat >> pyproject.toml <<EOF
[tool.meson-python]
limited-api = true
EOF
fi
