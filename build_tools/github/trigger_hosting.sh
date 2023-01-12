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

     if [[ "$PULL_REQUEST_NUMBER" == "null" ]]; then
          # The pull request is on the main (default) branch of the fork. The above API
          # call is unable to get the PR number associated with the commit:
          # https://docs.github.com/en/rest/commits/commits#list-pull-requests-associated-with-a-commit
          # We fallback to the search API here. The search API is not used everytime
          # because it has a lower rate limit.
          PULL_REQUEST_NUMBER=$(curl \
               -H "Accept: application/vnd.github+json" \
               -H "Authorization: token $GITHUB_TOKEN" \
               "https://api.github.com/search/issues?q=$COMMIT_SHA+repo:$GITHUB_REPOSITORY" 2>/dev/null \
               | jq '.items[0].number')
     fi

     BRANCH=pull/$PULL_REQUEST_NUMBER/head
else
     BRANCH=$HEAD_BRANCH
fi

# Circle CI REST API return HTTP response with 202 status code even when the POST requests fail.
# Hence we add some handling so that errors are reported on GitHub.
# For details see: https://circleci.com/docs/api/v2/index.html#operation/triggerPipeline
CIRCLE_CI_RESPONSE=$(curl --request POST \
     --url https://circleci.com/api/v2/project/gh/$GITHUB_REPOSITORY/pipeline \
     --header "Circle-Token: $CIRCLE_CI_TOKEN" \
     --header "content-type: application/json" \
     --header "x-attribution-actor-id: github_actions" \
     --header "x-attribution-login: github_actions" \
     --data \{\"branch\":\"$BRANCH\",\"parameters\":\{\"GITHUB_RUN_URL\":\"$GITHUB_RUN_URL\"\}\})

echo $CIRCLE_CI_RESPONSE

CIRCLE_CI_RESPONSE_MESSAGE=$(cat $CIRCLE_CI_RESPONSE | jq ".message")

if [[ "$CIRCLE_CI_RESPONSE_MESSAGE" == "null" ]]; then
  exit 0  # No message means there was no error.
fi

if [[ "$CIRCLE_CI_RESPONSE_MESSAGE" =~ .*"Not Found".* ]]; then
  echo "The endpoint has not been found on Circle CI REST API. Please check the request correctness."
fi

if [[ "$CIRCLE_CI_RESPONSE_MESSAGE" =~ .*"Permission denied".* ]]; then
  echo "Circle CI is blocking the start of the pipeline."
  echo "Please check for correct registration of tokens on Circle CI on: https://app.circleci.com/settings/user/tokens"
  echo "See: https://support.circleci.com/hc/en-us/articles/360050351292-How-to-Trigger-a-Workflow-via-CircleCI-API-v2"
fi

exit 1
