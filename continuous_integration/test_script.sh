#!/bin/bash

set -e

python --version
python -c "import numpy; print('numpy %s' % numpy.__version__)"
python -c "import scipy; print('scipy %s' % scipy.__version__)"
python setup.py build_ext --inplace

if [[ "$COVERAGE" == "true" ]]; then
    export WITH_COVERAGE="--with-coverage"
else
    export WITH_COVERAGE=""
fi
nosetests -s -v $WITH_COVERAGE sklearn
