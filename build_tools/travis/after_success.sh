#!/bin/bash
# This script is meant to be called by the "after_success" step defined in
# .travis.yml. See https://docs.travis-ci.com/ for more details.

# License: 3-clause BSD

set -e

if [[ "$COVERAGE" == "true" ]]; then
    # Need to run codecov from a git checkout, so we copy .coverage
    # from TEST_DIR where pytest has been run
    cp $TEST_DIR/.coverage $TRAVIS_BUILD_DIR
    cd $TRAVIS_BUILD_DIR
    # Ignore codecov failures as the codecov server is not
    # very reliable but we don't want travis to report a failure
    # in the github UI just because the coverage report failed to
    # be published.
    codecov || echo "codecov upload failed"
fi
