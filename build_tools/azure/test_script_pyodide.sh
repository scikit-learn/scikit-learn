#!/bin/bash

set -e

# We are using a pytest js wrapper script to run tests inside Pyodide. Maybe
# one day we can use a Pyodide venv instead but at the time of writing
# (2023-09-27) there is an issue with scipy.linalg in a Pyodide venv, see
# https://github.com/pyodide/pyodide/issues/3865 for more details.
node build_tools/azure/pytest-pyodide.js --pyargs sklearn --durations 20 --showlocals
