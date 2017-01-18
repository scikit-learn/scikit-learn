#!/bin/bash


# Add `master` branch to the update list.
# Otherwise CircleCI will give us a cached one.
FETCH_REFS="+master:master"

# Update PR refs for testing.
if [[ -n "${CIRCLE_PR_NUMBER}" ]]
then
    FETCH_REFS="${FETCH_REFS} +refs/pull/${CIRCLE_PR_NUMBER}/head:pr/${CIRCLE_PR_NUMBER}/head"
    FETCH_REFS="${FETCH_REFS} +refs/pull/${CIRCLE_PR_NUMBER}/merge:pr/${CIRCLE_PR_NUMBER}/merge"
fi

# Retrieve the refs.
git fetch -u origin ${FETCH_REFS}

# Checkout the PR merge ref.
if [[ -n "${CIRCLE_PR_NUMBER}" ]]
then
    git checkout -qf "pr/${CIRCLE_PR_NUMBER}/merge"
fi

# Check for merge conflicts.
if [[ -n "${CIRCLE_PR_NUMBER}" ]]
then
    git branch --merged | grep master > /dev/null
    git branch --merged | grep "pr/${CIRCLE_PR_NUMBER}/head" > /dev/null
fi
