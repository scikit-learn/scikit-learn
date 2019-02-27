##!/bin/bash

set -e

if [[ "$CHECK_PYTEST_SOFT_DEPENDENCY" == "true" ]]; then
    conda remove -y py pytest || pip uninstall -y py pytest
    if [[ "$COVERAGE" == "true" ]]; then
        # Need to append the coverage to the existing .coverage generated by
        # running the tests
        CMD="coverage run --append"
    else
        CMD="python"
    fi
    # .coverage from running the tests is in TEST_DIR
    cd $TEST_DIR
    $CMD -m sklearn.utils.tests.test_estimator_checks
    cd $OLDPWD
fi
