#!/bin/bash

set -e
pip install --upgrade pip
pip install pytest pytest-xdist

# Faster run of the source code tests
export USABLE_CPU_COUNT=`python -c "import joblib; print(joblib.cpu_count())"`
echo "Usable number of CPUs according to joblib: $USABLE_CPU_COUNT"
python -m pytest -n $USABLE_CPU_COUNT --pyargs sklearn

# Test that there are no links to system libraries
python -m threadpoolctl -i sklearn
