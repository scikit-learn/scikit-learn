#!/bin/bash

set -e
set -x

GITHUB_RUN_URL=https://nightly.link/$GITHUB_REPOSITORY/actions/runs/$RUN_ID

if [ "$EVENT" == pull_request ]
then
     PULL_REQUEST_NUMBER=$(curl \
          -H "Accept: application/vnd.github.v3+json" \
          -H "Authorization: token $GITHUB_TOKEN" \
          https://api.github.com/repos/$REPO_NAME/commits/$COMMIT_SHA/pulls 2>/dev/null \
          | jq '.[0].number')
     BRANCH=pull/$PULL_REQUEST_NUMBER/head
else
     BRANCH=$HEAD_BRANCH
fi

curl --request POST \
     --url https://circleci.com/api/v2/project/gh/$GITHUB_REPOSITORY/pipeline \
     --header "Circle-Token: $CIRCLE_CI_TOKEN" \
     --header "content-type: application/json" \
     --header "x-attribution-actor-id: github_actions" \
     --header "x-attribution-login: github_actions" \
     --data \{\"branch\":\"$BRANCH\",\"parameters\":\{\"GITHUB_RUN_URL\":\"$GITHUB_RUN_URL\"\}\}
