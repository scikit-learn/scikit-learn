#!/bin/bash

set -e

# Defines the show_installed_libraries and activate_environment functions.
source build_tools/shared.sh

activate_environment

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
    TEST_CMD="$TEST_CMD --cov-config='$COVERAGE_PROCESS_START' --cov sklearn --cov-report="
fi

if [[ -n "$CHECK_WARNINGS" ]]; then
    TEST_CMD="$TEST_CMD -Werror::DeprecationWarning -Werror::FutureWarning -Werror::sklearn.utils.fixes.VisibleDeprecationWarning"

    # Ignore pkg_resources deprecation warnings triggered by pyamg
    TEST_CMD="$TEST_CMD -W 'ignore:pkg_resources is deprecated as an API:DeprecationWarning'"
    TEST_CMD="$TEST_CMD -W 'ignore:Deprecated call to \`pkg_resources:DeprecationWarning'"

    # pytest-cov issue https://github.com/pytest-dev/pytest-cov/issues/557 not
    # fixed although it has been closed. https://github.com/pytest-dev/pytest-cov/pull/623
    # would probably fix it.
    TEST_CMD="$TEST_CMD -W 'ignore:The --rsyncdir command line argument and rsyncdirs config variable are deprecated.:DeprecationWarning'"

    # In some case, exceptions are raised (by bug) in tests, and captured by pytest,
    # but not raised again. This is for instance the case when Cython directives are
    # activated: IndexErrors (which aren't fatal) are raised on out-of-bound accesses.
    # In those cases, pytest instead raises pytest.PytestUnraisableExceptionWarnings,
    # which we must treat as errors on the CI.
    TEST_CMD="$TEST_CMD -Werror::pytest.PytestUnraisableExceptionWarning"

    # warnings has been fixed from dateutil main but not released yet, see
    # https://github.com/dateutil/dateutil/issues/1314
    TEST_CMD="$TEST_CMD -Wignore:datetime.datetime.utcfromtimestamp:DeprecationWarning"

    # Python 3.12 warnings from joblib fixed in master but not released yet,
    # see https://github.com/joblib/joblib/pull/1518
    TEST_CMD="$TEST_CMD -W 'ignore:ast.Num is deprecated:DeprecationWarning'"
    TEST_CMD="$TEST_CMD -W 'ignore:Attribute n is deprecated:DeprecationWarning'"

fi

if [[ "$PYTEST_XDIST_VERSION" != "none" && "$DISTRIB" != "conda-pip-scipy-dev" ]]; then
    XDIST_WORKERS=$(python -c "import joblib; print(joblib.cpu_count(only_physical_cores=True))")
    TEST_CMD="$TEST_CMD -n$XDIST_WORKERS"
fi

if [[ -n "$SELECTED_TESTS" ]]; then
    TEST_CMD="$TEST_CMD -k $SELECTED_TESTS"

    # Override to make selected tests run on all random seeds
    export SKLEARN_TESTS_GLOBAL_RANDOM_SEED="all"
fi

set -x
eval "$TEST_CMD --maxfail=10 --pyargs sklearn"
set +x
