#!/bin/bash

# The behavior of the script is controlled by environment variable
# defined in the docs.yml in the actions level folder of the project.

set -ex

DOC_REPO="scikit-learn.github.io"
BRANCH=$(basename "$GITHUB_REF")

if [ "$GITHUB_REF" =~ "main" ]
then
    DIR=dev
else
    # Strip off .X
    DIR="${BRANCH::-2}"
fi

MSG="Pushing the docs to $DIR/ for branch: $BRANCH, commit $GITHUB_SHA"

# check if it's a new branch
echo $DIR > .git/info/sparse-checkout

if ! git show HEAD:$DIR >/dev/null
then
	# directory does not exist. Need to make it so sparse checkout works
	mkdir $DIR
	touch $DIR/index.html
	git add $DIR
fi
git checkout main
git reset --hard origin/main
if [ -d $DIR ]
then
	git rm -rf $DIR/ && rm -rf $DIR/
fi

cp -R ../docs $DIR

ls -l ../docs
ls -l $DIR

git config user.email "olivier.grisel+sklearn-ci@gmail.com"
git config user.name "sklearn-ci"
git config push.default matching

git add -f $DIR/
git commit -m "$MSG" $DIR
git push

echo $MSG
