#!/bin/bash

set -e

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

python -c "import sklearn; sklearn.show_versions()"
python -m threadpoolctl -i sklearn

if ! command -v conda &> /dev/null
then
    pip list
else
    # conda list provides more info than pip list (when available)
    conda list
fi

TEST_CMD="python -m pytest --showlocals --durations=20 --junitxml=$JUNITXML"

if [[ "$COVERAGE" == "true" ]]; then
    # Note: --cov-report= is used to disable to long text output report in the
    # CI logs. The coverage data is consolidated by codecov to get an online
    # web report across all the platforms so there is no need for this text
    # report that otherwise hides the test failures and forces long scrolls in
    # the CI logs.
    export COVERAGE_PROCESS_START="$BUILD_SOURCESDIRECTORY/.coveragerc"
    TEST_CMD="$TEST_CMD --cov-config=$COVERAGE_PROCESS_START --cov sklearn --cov-report="
fi

if [[ -n "$CHECK_WARNINGS" ]]; then
    # numpy's 1.19.0's tostring() deprecation is ignored until scipy and joblib removes its usage
    TEST_CMD="$TEST_CMD -Werror::DeprecationWarning -Werror::FutureWarning -Wignore:tostring:DeprecationWarning"
fi

if [[ "$PYTEST_XDIST_VERSION" != "none" ]]; then
    TEST_CMD="$TEST_CMD -n2"
fi

set -x
$TEST_CMD --pyargs sklearn
set +x
