#!/bin/bash

set -e
set -x

pytest --pyargs sklearn

# Test that there are no links to system libraries
python -m threadpoolctl -i sklearn
