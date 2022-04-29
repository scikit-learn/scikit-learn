#!/bin/bash

set -e
set -x

GITHUB_RUN_URL=https://nightly.link/$GITHUB_REPOSITORY/actions/runs/$RUN_ID
BRANCH=$( wget $GITHUB_RUN_URL/branch.zip && \
          unzip branch.zip > /dev/null && \
          cat branch.txt )

curl --request POST \
     --url https://circleci.com/api/v2/project/gh/$GITHUB_REPOSITORY/pipeline \
     --header "Circle-Token: $CCI_TOKEN" \
     --header "content-type: application/json" \
     --header "x-attribution-actor-id: github_actions" \
     --header "x-attribution-login: github_actions" \
     --data \{\"branch\":\"$BRANCH\",\"parameters\":\{\"GITHUB_RUN_URL\":$GITHUB_RUN_URL\"\"\}\}
