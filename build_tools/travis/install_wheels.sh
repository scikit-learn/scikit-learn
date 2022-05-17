#!/bin/bash

python -m pip install cibuildwheel || travis_terminate $?
python -m cibuildwheel --output-dir wheelhouse || travis_terminate $?
