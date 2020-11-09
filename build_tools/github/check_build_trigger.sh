#!/bin/bash

set -e
set -x

COMMIT_MSG=$(git log -1 --pretty=%B)

# By default pull requests use refs/pull/PULL_ID/merge as the source branch
# which has a "Merge ID into ID" as a commit message. The latest commit
# message is the second to last commit
if [ "$GITHUB_EVENT_NAME" == "pull_request" ]; then
  COMMIT_ID=$(echo $MERGE_MSG | awk '{print $2}')
  COMMIT_MSG=$(git log $COMMIT_ID -1 --pretty=%B)
fi

# The commit marker "[cd build]" will trigger the build when required
if [[ "$GITHUB_EVENT_NAME" == schedule ||
      "$COMMIT_MSG" =~ \[cd\ build\] ]]; then
    echo "::set-output name=build::true"
fi
