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

pyodide venv pyodide-venv
source pyodide-venv/bin/activate

pip install dist/*.whl
pip list
