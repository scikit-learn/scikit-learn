#!/bin/bash

# License: 3-clause BSD

# This script is meant to be called by the "after_success" step
# defined in ".travis.yml". In particular, we upload the wheels
# for the continuous deployment jobs and the code coverage for
# the continuous integration jobs.

set -e
set -x

if [[ $BUILD_WHEEL == true ]]; then
    # TODO: Upload the wheels to the Anaconda repository
else
    if [[ "$COVERAGE" == "true" ]]; then
        # Need to run codecov from a git checkout, so we copy .coverage
        # from TEST_DIR where pytest has been run
        cp $TEST_DIR/.coverage $TRAVIS_BUILD_DIR

        # Ignore codecov failures as the codecov server is not
        # very reliable but we don't want travis to report a failure
        # in the github UI just because the coverage report failed to
        # be published.
        codecov --root $TRAVIS_BUILD_DIR || echo "codecov upload failed"
    fi
fi
