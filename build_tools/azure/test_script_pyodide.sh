#!/bin/bash

set -e

source pyodide-venv/bin/activate

pip list

# Need to be outside of the git clone otherwise finds non build sklearn folder
cd /tmp

# TODO for now only testing sklearn import to make sure the wheel is not badly
# broken. When Pyodide 0.24 is released we should run the full test suite and
# xfail tests that fail due to Pyodide limitations
python -c 'import sklearn'
