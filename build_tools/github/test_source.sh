#!/bin/bash

set -e
set -x

pytest --pyargs sklearn

# Check whether the source distribution will render correctly
twine check dist/*.tar.gz
