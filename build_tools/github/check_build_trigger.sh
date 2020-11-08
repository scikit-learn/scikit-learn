#!/bin/bash

set -e
set -x

COMMIT_ID=$(echo $GITHUB_SHA | awk '{print $2}')
COMMIT_MESSAGE=$(git log $COMMIT_ID -1 --pretty=%B)

# The commit marker "[cd build]" will trigger the build when required
if [[ "$GITHUB_EVENT_NAME" == schedule ||
      "$COMMIT_MSG" =~ \[cd\ build\] ]]; then
    echo "::set-output name=build::true"
fi
