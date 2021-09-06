#!/bin/bash

set -e
pip install --upgrade pip
pip install pytest pytest-xdist

# Faster run of the source code tests
python -m pytest -n $CPU_COUNT --pyargs sklearn

# Test that there are no links to system libraries
python -m threadpoolctl -i sklearn
