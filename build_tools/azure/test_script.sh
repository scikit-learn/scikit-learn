#!/bin/bash

set -e

if [[ "$DISTRIB" == "conda" ]]; then
    source activate $VIRTUALENV
elif [[ "$DISTRIB" == "ubuntu" ]]; then
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

TEST_CMD="python -m pytest --showlocals --durations=20 --junitxml=$JUNITXML --pyargs"

if [[ "$COVERAGE" == "true" ]]; then
    TEST_CMD="$TEST_CMD --cov sklearn"
fi

if [[ -n "$CHECK_WARNINGS" ]]; then
    TEST_CMD="$TEST_CMD -Werror::DeprecationWarning -Werror::FutureWarning"
fi

mkdir -p $TEST_DIR
cp setup.cfg $TEST_DIR
cd $TEST_DIR

# Collects tests 3 times and confirms that they match
# This is to make sure the tests will continue to work with pytest-xdist.
if [[ "$PYTHON_VERSION" == "*" ]]; then
    pytest --collect-only -q --pyargs $1 | sed '$d' > collect_1.txt
    pytest --collect-only -q --pyargs $1 | sed '$d' > collect_2.txt
    diff collect_2.txt collect_1.txt

    pytest --collect-only -q --pyargs $1 | sed '$d' > collect_3.txt
    diff collect_3.txt collect_1.txt
    diff collect_3.txt collect_2.txt
fi

set -x
$TEST_CMD sklearn
set +x
