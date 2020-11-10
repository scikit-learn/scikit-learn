#!/bin/bash

set -e
set -x

python -m pip install cibuildwheel

# Extend the waiting time to avoid that the tests stall
travis_wait 30 python -m cibuildwheel --output-dir wheelhouse
