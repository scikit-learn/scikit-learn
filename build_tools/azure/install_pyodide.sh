#!/bin/bash

set -e

git clone https://github.com/emscripten-core/emsdk.git
cd emsdk
./emsdk install $EMSCRIPTEN_VERSION
./emsdk activate $EMSCRIPTEN_VERSION
source emsdk_env.sh
cd -

pip install pyodide-build==$PYODIDE_VERSION pyodide-cli

# TODO: temporary fix for pywasmcross.py taken from https://github.com/pyodide/pyodide/pull/4136
wget https://github.com/lesteve/pyodide/raw/pywasmcross-include/pyodide-build/pyodide_build/pywasmcross.py
PYODIDE_BUILD_DIR=$(python -c 'import os; import pyodide_build; print(os.path.dirname(pyodide_build.__file__))')
echo ${PYODIDE_BUILD_DIR}
ls ${PYODIDE_BUILD_DIR}
diff -u pywasmcross.py ${PYODIDE_BUILD_DIR}/pywasmcross.py
cp -f pywasmcross.py ${PYODIDE_BUILD_DIR}/pywasmcross.py
diff -u pywasmcross.py ${PYODIDE_BUILD_DIR}/pywasmcross.py

pyodide build

ls -ltrh dist

pyodide venv pyodide-venv
source pyodide-venv/bin/activate

pip install dist/*.whl
pip list
