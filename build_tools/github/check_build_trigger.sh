#!/bin/bash

set -e
set -x

COMMIT_MSG=$(git log --no-merges -1 --oneline)

ON_PUSH = $GITHUB_EVENT_NAME == "push"
ON_SCHEDULE = $GITHUB_EVENT_NAME == "schedule"

# The commit marker "[cd build]" will trigger the build when required
ON_COMMIT_MARKER = $COMMIT_MSG =~ "[cd build]"

if [ $ON_PUSH | $ON_SCHEDULE | $ON_COMMIT_MARKER ]; then
    echo "::set-output name=build::true"
fi
