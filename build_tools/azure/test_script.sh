#!/bin/bash

set -e

# defines the show_installed_libraries function
source build_tools/shared.sh

if [[ "$DISTRIB" =~ ^conda.* ]]; then
    source activate $VIRTUALENV
elif [[ "$DISTRIB" == "ubuntu" ]] || [[ "$DISTRIB" == "debian-32" ]]; then
    source $VIRTUALENV/bin/activate
fi

if [[ "$BUILD_WITH_ICC" == "true" ]]; then
    source /opt/intel/oneapi/setvars.sh
fi

mkdir -p $TEST_DIR
cp setup.cfg $TEST_DIR
cd $TEST_DIR

# sklearn.show_versions() now strictly requires Python 3.8+
# TODO: remove this once PyPy 3.8+ is released on conda-forge
PYTHON_IMPL=`python -c "import platform; print(platform.python_implementation())"`
if [[ "$PYTHON_IMPL" != "PyPy" ]]; then
    python -c "import sklearn; sklearn.show_versions()"
fi

show_installed_libraries

TEST_CMD="python -m pytest --showlocals --durations=20 --junitxml=$JUNITXML"

if [[ "$COVERAGE" == "true" ]]; then
    # Note: --cov-report= is used to disable to long text output report in the
    # CI logs. The coverage data is consolidated by codecov to get an online
    # web report across all the platforms so there is no need for this text
    # report that otherwise hides the test failures and forces long scrolls in
    # the CI logs.
    export COVERAGE_PROCESS_START="$BUILD_SOURCESDIRECTORY/.coveragerc"
    TEST_CMD="$TEST_CMD --cov-config='$COVERAGE_PROCESS_START' --cov sklearn --cov-report="
fi

if [[ -n "$CHECK_WARNINGS" ]]; then
    # numpy's 1.19.0's tostring() deprecation is ignored until scipy and joblib removes its usage
    TEST_CMD="$TEST_CMD -Werror::DeprecationWarning -Werror::FutureWarning -Wignore:tostring:DeprecationWarning"

    # Python 3.10 deprecates disutils and is imported by numpy interally during import time
    TEST_CMD="$TEST_CMD -Wignore:The\ distutils:DeprecationWarning"
fi

if [[ "$PYTEST_XDIST_VERSION" != "none" ]]; then
    TEST_CMD="$TEST_CMD -n$CPU_COUNT"
fi

if [[ "$SHOW_SHORT_SUMMARY" == "true" ]]; then
    TEST_CMD="$TEST_CMD -ra"
fi

set -x
eval "$TEST_CMD --pyargs sklearn"
set +x
