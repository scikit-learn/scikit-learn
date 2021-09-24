#!/bin/bash

# The behavior of the script is controlled by environment variable
# defined in the docs.yml in the actions level folder of the project.

set -ex

cd scikit-learn.github.io

DOC_REPO="scikit-learn.github.io"

if [[ -z "$GITHUB_BASE_REF" ]]
then
	REF="$GITHUB_REF"
else
	REF="$GITHUB_BASE_REF"
fi

if [[ "$REF" =~ "main" ]]
then
    DIR=dev
else
    # Strip off .X
    DIR="${REF::-2}"
fi

MSG="Pushing the docs to $DIR/ for branch: $REF, commit $GITHUB_SHA"

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

ls -l $GITHUB_WORKSPACE/docs

cp -R $GITHUB_WORKSPACE/docs $DIR

git config user.email "olivier.grisel+sklearn-ci@gmail.com"
git config user.name "sklearn-ci"
git config push.default matching

git add -f $DIR/
git commit -m "$MSG" $DIR
# git push

echo $MSG
