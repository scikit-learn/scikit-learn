#!/bin/bash

set -e
set -x

cp $CONFTEST_PATH $CONFTEST_NAME

pytest --pyargs sklearn

# Test that there are no links to system libraries
python -m threadpoolctl -i sklearn
