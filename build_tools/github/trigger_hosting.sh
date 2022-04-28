#!/bin/bash

set -e
set -x

curl --request POST \
     --url https://circleci.com/api/v2/project/gh/$GITHUB_REPOSITORY/continue \
     --header "Circle-Token: $CCI_TOKEN" \
     --header "content-type: application/json" \
     --header "x-attribution-actor-id: github_actions" \
     --header "x-attribution-login: github_actions" \
     --data \{\"parameters\":\{\"GITHUB_RUN_URL\":\"$GITHUB_RUN_URL\"\}\}
