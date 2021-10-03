#!/bin/bash

set -e
set -x

# COMMIT_MSG=$(git log --no-merges -1 --oneline)
COMMIT_MSG="wow"

# The commit marker "[cd build]" will trigger the build when required
if [[ "$GITHUB_EVENT_NAME" == schedule ||
      "$COMMIT_MSG" =~ \[cd\ build\] ||
      "$COMMIT_MSG" =~ \[cd\ build\ gh\] ]]; then
    echo "::set-output name=build::true"
fi
