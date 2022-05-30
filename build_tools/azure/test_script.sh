#!/bin/bash

set -e

# defines the show_installed_libraries function
source build_tools/shared.sh

if [[ "$DISTRIB" =~ ^conda.* ]]; then
    source activate $VIRTUALENV
elif [[ "$DISTRIB" == "ubuntu" || "$DISTRIB" == "debian-32" || "$DISTRIB" == "pip-nogil" ]]; then
    source $VIRTUALENV/bin/activate
elif [[ "$DISTRIB" == "pip-windows" ]]; then
    source $VIRTUALENV/Scripts/activate
fi

if [[ "$BUILD_WITH_ICC" == "true" ]]; then
    source /opt/intel/oneapi/setvars.sh
fi

if [[ "$BUILD_REASON" == "Schedule" ]]; then
    # Enable global random seed randomization to discover seed-sensitive tests
    # only on nightly builds.
    # https://scikit-learn.org/stable/computing/parallelism.html#environment-variables
    export SKLEARN_TESTS_GLOBAL_RANDOM_SEED="any"

    # Enable global dtype fixture for all nightly builds to discover
    # numerical-sensitive tests.
    # https://scikit-learn.org/stable/computing/parallelism.html#environment-variables
    export SKLEARN_RUN_FLOAT32_TESTS=1
fi

mkdir -p $TEST_DIR
cp setup.cfg $TEST_DIR
cd $TEST_DIR

python -c "import joblib; print(f'Number of cores: {joblib.cpu_count()}')"
python -c "import sklearn; sklearn.show_versions()"

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
    TEST_CMD="$TEST_CMD -Werror::DeprecationWarning -Werror::FutureWarning"

    # numpy's 1.19.0's tostring() deprecation is ignored until scipy and joblib
    # removes its usage
    TEST_CMD="$TEST_CMD -Wignore:tostring:DeprecationWarning"

    # Python 3.10 deprecates distutils, which is imported by numpy internally
    TEST_CMD="$TEST_CMD -Wignore:The\ distutils:DeprecationWarning"

    # Ignore distutils deprecation warning, used by joblib internally
    TEST_CMD="$TEST_CMD -Wignore:distutils\ Version\ classes\ are\ deprecated:DeprecationWarning"
fi

if [[ "$PYTEST_XDIST_VERSION" != "none" ]]; then
    TEST_CMD="$TEST_CMD -n$CPU_COUNT"
fi

if [[ "$SHOW_SHORT_SUMMARY" == "true" ]]; then
    TEST_CMD="$TEST_CMD -ra"
fi

if [[ -n "$SELECTED_TESTS" ]]; then
    TEST_CMD="$TEST_CMD -k $SELECTED_TESTS"

    # Override to make selected tests run on all random seeds
    export SKLEARN_TESTS_GLOBAL_RANDOM_SEED="all"
fi

set -x
eval "$TEST_CMD --pyargs sklearn"
set +x
