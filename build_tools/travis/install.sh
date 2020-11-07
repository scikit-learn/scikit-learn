#!/bin/bash

# This script is meant to be called by the "install" step
# defined in the ".travis.yml" file. In particular, it is
# important that we call to the right installation script
# for the job being executed.

set -e
set -x

# We cannot use the "TRAVIS_COMMIT_MESSAGE" environment
# variable because for PRs it contains the merge commit 
COMMIT_MSG=$(git log --no-merges -1 --oneline)

if [[ "$COMMIT_MSG" =~ \[cd\ build\] ]]; then
    source build_tools/travis/install_wheels.sh
else
    source build_tools/travis/install_master.sh
fi
