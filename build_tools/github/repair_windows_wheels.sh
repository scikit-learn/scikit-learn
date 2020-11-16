#!/bin/bash

set -e
set -x

WHEEL=$1
DEST_DIR=$2

# By default, the Windows wheels are not repaired.
# In this case, we want to vendor the vcomp140.dll
wheel unpack "$WHEEL"
python build_tools/github/vendor_vcomp140.py "$WHEEL_DIRNAME"
wheel pack "$WHEEL_DIRNAME" -d "$DEST_DIR"
rm -rf "$WHEEL_DIRNAME"
