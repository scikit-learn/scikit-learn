#!/bin/bash

set -e
set -x

WHEEL=$1
DEST_DIR=$2
BITNESS=$3

# By default, the Windows wheels are not repaired. In this case,
# it is neccesary to vendor vcomp140.dll and msvcp140.dll files
wheel unpack "$WHEEL"
python build_tools/github/vendor.py "$WHEEL_DIRNAME" "$BITNESS"
wheel pack "$WHEEL_DIRNAME" -d "$DEST_DIR"
rm -rf "$WHEEL_DIRNAME"
