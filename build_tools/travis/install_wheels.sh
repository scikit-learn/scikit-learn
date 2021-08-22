#!/bin/bash

set -e

python -m pip install cibuildwheel
python -m cibuildwheel --output-dir wheelhouse
