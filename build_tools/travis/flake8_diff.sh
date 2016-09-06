#!/bin/bash

# This script is used in Travis to check that PRs do not add obvious
# flake8 violations. It relies on two things:
#   - find common ancestor between branch and
#     scikit-learn/scikit-learn remote
#   - run flake8 --diff on the diff between the branch and the common
#     ancestor
#
# Additional features:
#   - the line numbers in Travis match the local branch on the PR
#     author machine.
#   - ./build_tools/travis/flake8_diff.sh can be run locally for quick
#     turn-around

set -e
# pipefail is necessary to propagate exit codes
set -o pipefail

PROJECT=scikit-learn/scikit-learn
PROJECT_URL=https://github.com/$PROJECT.git

echo "Remotes:"
git remote --verbose

# Find the remote with the project name (upstream in most cases)
REMOTE=$(git remote -v | grep $PROJECT | cut -f1 | head -1 || echo '')

# Add a temporary remote if needed. For example this is necessary when
# Travis is configured to run in a fork. In this case 'origin' is the
# fork and not the reference repo we want to diff against.
if [[ -z "$REMOTE" ]]; then
    TMP_REMOTE=tmp_reference_upstream
    REMOTE=$TMP_REMOTE
    git remote add $REMOTE $PROJECT_URL
fi

if [[ "$TRAVIS" == "true" ]]; then
    if [[ "$TRAVIS_PULL_REQUEST" == "false" ]]
    then
        # Travis does the git clone with a limited depth (50 at the time of
        # writing). This may not be enough to find the common ancestor with
        # $REMOTE/master so we unshallow the git checkout
        git fetch --unshallow || echo "Unshallowing the git checkout failed"
    else
        # We want to fetch the code as it is in the PR branch and not
        # the result of the merge into master. This way line numbers
        # reported by Travis will match with the local code.
        BRANCH_NAME=travis_pr_$TRAVIS_PULL_REQUEST
        git fetch $REMOTE pull/$TRAVIS_PULL_REQUEST/head:$BRANCH_NAME
        git checkout $BRANCH_NAME
    fi
fi


echo -e '\nLast 2 commits:'
echo '--------------------------------------------------------------------------------'
git log -2 --pretty=short

git fetch $REMOTE master
REMOTE_MASTER_REF="$REMOTE/master"

# Find common ancestor between HEAD and remotes/$REMOTE/master
COMMIT=$(git merge-base @ $REMOTE_MASTER_REF) || \
    echo "No common ancestor found for $(git show @ -q) and $(git show $REMOTE_MASTER_REF -q)"

if [[ -n "$TMP_REMOTE" ]]; then
    git remote remove $TMP_REMOTE
fi

if [ -z "$COMMIT" ]; then
    exit 1
fi

echo -e "\nCommon ancestor between HEAD and $REMOTE_MASTER_REF is:"
echo '--------------------------------------------------------------------------------'
git show --no-patch $COMMIT

echo -e '\nRunning flake8 on the diff in the range'\
     "$(git rev-parse --short $COMMIT)..$(git rev-parse --short @)" \
     "($(git rev-list $COMMIT.. | wc -l) commit(s)):"
echo '--------------------------------------------------------------------------------'

# We ignore files from sklearn/externals. Unfortunately there is no
# way to do it with flake8 directly (the --exclude does not seem to
# work with --diff). We could use the exclude magic in the git pathspec
# ':!sklearn/externals' but it is only available on git 1.9 and Travis
# uses git 1.8.
# We need the following command to exit with 0 hence the echo in case
# there is no match
MODIFIED_FILES=$(git diff --name-only $COMMIT | grep -v 'sklearn/externals' || echo "no_match")

if [[ "$MODIFIED_FILES" == "no_match" ]]; then
    echo "No file outside sklearn/externals has been modified"
else
    # Conservative approach: diff without context so that code that
    # was not changed does not create failures
    git diff --unified=0 $COMMIT -- $MODIFIED_FILES | flake8 --diff --show-source
fi
echo -e "No problem detected by flake8\n"
