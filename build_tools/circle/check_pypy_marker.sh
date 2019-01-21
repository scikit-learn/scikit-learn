#!/usr/bin/env bash
set -x
set -e

# Check for pypy marker in the commit message and
# trigger the job if the marker is present.

commit_msg=$(git log --format=%B -n 1 $CIRCLE_SHA1)

if [[ "$commit_msg" =~ "[pypy]" ]]
    then
        curl --user ${CIRCLE_API_USER_TOKEN} \
            --data build_parameters[CIRCLE_JOB]=pypy3 \
            https://circleci.com/api/v1.1/project/github/$CIRCLE_PROJECT_USERNAME/$CIRCLE_PROJECT_REPONAME/tree/$CIRCLE_BRANCH
fi
