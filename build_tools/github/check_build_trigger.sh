#!/bin/bash

set -e
set -x

COMMIT_MSG=$(git log --no-merges -1 --oneline)
BRANCH_NAME=$(echo "${GITHUB_BRANCH##*/}")

# The commit marker "[cd build]" will trigger the build when required
if [[ "$GITHUB_EVENT_NAME" == push && "$BRANCH_NAME" != master ||
      "$GITHUB_EVENT_NAME" == schedule ||
      "$COMMIT_MSG" =~ \[cd\ build\] ]]; then
    echo "::set-output name=build::true"
fi
