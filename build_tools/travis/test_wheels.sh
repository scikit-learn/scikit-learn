#!/bin/bash

set -e

# Faster run of the source code tests
pytest -n $CPU_COUNT --pyargs sklearn

# Test that there are no links to system libraries
python -m threadpoolctl -i sklearn
