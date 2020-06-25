#!/bin/bash

set -e

if [[ "$DISTRIB" =~ ^conda.* ]]; then
    source activate $VIRTUALENV
elif [[ "$DISTRIB" == "ubuntu" ]] || [[ "$DISTRIB" == "ubuntu-32" ]]; then
    source $VIRTUALENV/bin/activate
fi

python --version
python -c "import numpy; print('numpy %s' % numpy.__version__)"
python -c "import scipy; print('scipy %s' % scipy.__version__)"
python -c "\
try:
    import pandas
    print('pandas %s' % pandas.__version__)
except ImportError:
    print('pandas not installed')
"
python -c "import multiprocessing as mp; print('%d CPUs' % mp.cpu_count())"
pip list

TEST_CMD="python -m pytest --showlocals --durations=20 --junitxml=$JUNITXML"

if [[ "$COVERAGE" == "true" ]]; then
    export COVERAGE_PROCESS_START="$BUILD_SOURCESDIRECTORY/.coveragerc"
    TEST_CMD="$TEST_CMD --cov-config=$COVERAGE_PROCESS_START --cov sklearn"
fi

if [[ -n "$CHECK_WARNINGS" ]]; then
    TEST_CMD="$TEST_CMD -Werror::DeprecationWarning -Werror::FutureWarning"
    # numpy's 1.19.0's tostring() deprecation is ignored until scipy and joblib removes its usage
    TEST_CMD="$TEST_CMD -Wignore:tostring:DeprecationWarning"
    # numpy's 1.20.0's np.bool, np.int, np.float and np.object deprecations are ignored until pandas removes its usage
    TEST_CMD="$TEST_CMD -Wignore:\`np.*\`:DeprecationWarning:pandas"
fi

if [[ "$PYTEST_XDIST" == "true" ]]; then
    TEST_CMD="$TEST_CMD -n2"
fi

mkdir -p $TEST_DIR
cp setup.cfg $TEST_DIR
cd $TEST_DIR

set -x
$TEST_CMD --pyargs sklearn
set +x
