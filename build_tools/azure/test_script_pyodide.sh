#!/bin/bash

set -e

npm install pyodide@$PYODIDE_VERSION

node build_tools/azure/pytest-pyodide.js --pyargs sklearn
