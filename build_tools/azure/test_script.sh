#!/bin/bash

set -e

UNAMESTR=`uname`

python --version
python -c "import numpy; print('numpy %s' % numpy.__version__)"
python -c "import scipy; print('scipy %s' % scipy.__version__)"
python -c "import multiprocessing as mp; print('%d CPUs' % mp.cpu_count())"
python -c "\
try:
    import pandas
    print('pandas %s' % pandas.__version__)
except ImportError:
    pass
"

run_tests() {
    TEST_CMD="pytest --showlocals --durations=20 --pyargs"

    mkdir -p $TEST_DIR
    cp setup.cfg $TEST_DIR
    cd $TEST_DIR

    export SKLEARN_SKIP_NETWORK_TESTS=1

    set -x
    $TEST_CMD sklearn/preprocessing/tests/
}


if [[ "$UNAMESTR" == "Linux" ]]; then
    source $VIRTUALENV_DIR/bin/activate
fi

echo $TEST_DIR
# run_tests


