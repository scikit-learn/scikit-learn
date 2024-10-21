#!/bin/bash

set -e

git clone https://github.com/emscripten-core/emsdk.git
cd emsdk
./emsdk install $EMSCRIPTEN_VERSION
./emsdk activate $EMSCRIPTEN_VERSION
source emsdk_env.sh
cd -

pip install pyodide-build==$PYODIDE_VERSION pyodide-cli

pyodide build

ls -ltrh dist

# The Pyodide js library is needed by build_tools/azure/test_script_pyodide.sh
# to run tests inside Pyodide
npm install pyodide@$PYODIDE_VERSION
