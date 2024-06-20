#!/bin/bash

set -e
set -x

COMMIT_MSG=$(git log --no-merges -1 --oneline)

# The commit marker "[cd build]" or "[cd build gh]" will trigger the build when required
if [[ "$GITHUB_EVENT_NAME" == schedule ||
      "$GITHUB_EVENT_NAME" == workflow_dispatch ||
      "$COMMIT_MSG" =~ \[cd\ build\] ||
      "$COMMIT_MSG" =~ \[cd\ build\ gh\] ]]; then
    echo "build=true" >> $GITHUB_OUTPUT
fi
