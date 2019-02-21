#!/bin/bash

set -e

if [[ "$COVERAGE" == "true" ]]; then
    if [[ "$DISTRIB" == "conda" ]]; then
        export PATH=$HOME/miniconda3/bin:$PATH
        source activate $VIRTUALENV
    elif [[ "$DISTRIB" == "ubuntu" ]]; then
        source $VIRTUALENV/bin/activate
    fi
    # Need to run codecov from a git checkout, so we copy .coverage
    # from TEST_DIR where pytest has been run
    cp $TEST_DIR/.coverage $TRAVIS_BUILD_DIR

    # Ignore codecov failures as the codecov server is not
    # very reliable but we don't want travis to report a failure
    # in the github UI just because the coverage report failed to
    # be published.
    codecov --root $TRAVIS_BUILD_DIR || echo "codecov upload failed"
fi
