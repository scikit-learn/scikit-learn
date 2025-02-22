#!/bin/bash

set -e

# Defines the show_installed_libraries and activate_environment functions.
source build_tools/shared.sh

activate_environment

if [[ "$BUILD_REASON" == "Schedule" ]]; then
    # Enable global random seed randomization to discover seed-sensitive tests
    # only on nightly builds.
    # https://scikit-learn.org/stable/computing/parallelism.html#environment-variables
    export SKLEARN_TESTS_GLOBAL_RANDOM_SEED=$(($RANDOM % 100))
    echo "To reproduce this test run, set the following environment variable:"
    echo "    SKLEARN_TESTS_GLOBAL_RANDOM_SEED=$SKLEARN_TESTS_GLOBAL_RANDOM_SEED",
    echo "See: https://scikit-learn.org/dev/computing/parallelism.html#sklearn-tests-global-random-seed"

    # Enable global dtype fixture for all nightly builds to discover
    # numerical-sensitive tests.
    # https://scikit-learn.org/stable/computing/parallelism.html#environment-variables
    export SKLEARN_RUN_FLOAT32_TESTS=1
fi

COMMIT_MESSAGE=$(python build_tools/azure/get_commit_message.py --only-show-message)

if [[ "$COMMIT_MESSAGE" =~ \[float32\] ]]; then
    echo "float32 tests will be run due to commit message"
    export SKLEARN_RUN_FLOAT32_TESTS=1
fi

mkdir -p $TEST_DIR
cp setup.cfg $TEST_DIR
cd $TEST_DIR

python -c "import joblib; print(f'Number of cores (physical): \
{joblib.cpu_count()} ({joblib.cpu_count(only_physical_cores=True)})')"
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

    # Use sys.monitoring to make coverage faster for Python >= 3.12
    HAS_SYSMON=$(python -c 'import sys; print(sys.version_info >= (3, 12))')
    if [[ "$HAS_SYSMON" == "True" ]]; then
        export COVERAGE_CORE=sysmon
    fi
    TEST_CMD="$TEST_CMD --cov-config='$COVERAGE_PROCESS_START' --cov sklearn --cov-report="
fi

if [[ "$PYTEST_XDIST_VERSION" != "none" ]]; then
    XDIST_WORKERS=$(python -c "import joblib; print(joblib.cpu_count(only_physical_cores=True))")
    TEST_CMD="$TEST_CMD -n$XDIST_WORKERS"
fi

if [[ -n "$SELECTED_TESTS" ]]; then
    TEST_CMD="$TEST_CMD -k $SELECTED_TESTS"

    # Override to make selected tests run on all random seeds
    export SKLEARN_TESTS_GLOBAL_RANDOM_SEED="all"
fi

if which lscpu ; then
    lscpu
else
    echo "Could not inspect CPU architecture."
fi

TEST_CMD="$TEST_CMD --pyargs sklearn"

set -x
eval "$TEST_CMD"
set +x
