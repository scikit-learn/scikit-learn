#!/bin/bash

pip install --upgrade pip || travis_terminate $?
pip install pytest pytest-xdist || travis_terminate $?

# Faster run of the source code tests
export USABLE_CPU_COUNT=`python -c "import joblib; print(joblib.cpu_count())"`
echo "Usable number of CPUs according to joblib: $USABLE_CPU_COUNT"
python -m pytest -n $USABLE_CPU_COUNT --pyargs sklearn || travis_terminate $?

# Test that there are no links to system libraries
python -m threadpoolctl -i sklearn || travis_terminate $?
