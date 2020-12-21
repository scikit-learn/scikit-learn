#!/bin/bash

set -e

python --version
python -c "import numpy; print(f'numpy {numpy.__version__}')"
python -c "import scipy; print(f'scipy {scipy.__version__}')"
python -c "\
try:
    import pandas
    print(f'pandas {pandas.__version__}')
except ImportError:
    pass
"
python -c "import joblib; print(f'{joblib.cpu_count()} CPUs')"
python -c "import platform; print(f'{platform.machine()}')"

TEST_CMD="pytest --showlocals --durations=20 --pyargs"

# Run the tests on the installed version
mkdir -p $TEST_DIR

# Copy "setup.cfg" for the test settings
cp setup.cfg $TEST_DIR
cd $TEST_DIR

if [[ $TRAVIS_CPU_ARCH == arm64 ]]; then
    # Faster run of the source code tests
    TEST_CMD="$TEST_CMD -n $CPU_COUNT"

    # Remove the option to test the docstring
    sed -i -e 's/--doctest-modules//g' setup.cfg
fi

if [[ -n $CHECK_WARNINGS ]]; then
    TEST_CMD="$TEST_CMD -Werror::DeprecationWarning -Werror::FutureWarning"
fi

$TEST_CMD sklearn
