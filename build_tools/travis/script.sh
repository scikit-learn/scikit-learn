#!/bin/bash

# This script is meant to be called by the "script" step defined
# in the ".travis.yml" file. While this step is forbidden by the
# continuous deployment jobs, we execute the testing scripts for
# the continuous integration jobs.

set -e
set -x

if ! [[ "$TRAVIS_COMMIT_MESSAGE" =~ \[cd\ build\] ]]; then
    bash build_tools/travis/test_script.sh || travis_terminate 1
    bash build_tools/travis/test_docs.sh || travis_terminate 1
    bash build_tools/travis/test_pytest_soft_dependency.sh || travis_terminate 1
fi
