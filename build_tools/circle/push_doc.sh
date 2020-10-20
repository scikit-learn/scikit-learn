#!/bin/bash
# This script is meant to be called in the "deploy" step defined in 
# circle.yml. See https://circleci.com/docs/ for more details.
# The behavior of the script is controlled by environment variable defined
# in the circle.yml in the top level folder of the project.

set -ex

if [ -z $CIRCLE_PROJECT_USERNAME ];
then USERNAME="sklearn-ci";
else USERNAME=$CIRCLE_PROJECT_USERNAME;
fi

DOC_REPO="scikit-learn.github.io"
GENERATED_DOC_DIR=$1

if [[ -z "$GENERATED_DOC_DIR" ]]; then
    echo "Need to pass directory of the generated doc as argument"
    echo "Usage: $0 <generated_doc_dir>"
    exit 1
fi

# Absolute path needed because we use cd further down in this script
GENERATED_DOC_DIR=$(readlink -f $GENERATED_DOC_DIR)

if [ "$CIRCLE_BRANCH" = "master" ]
then
    dir=dev
else
    # Strip off .X
    dir="${CIRCLE_BRANCH::-2}"
fi

MSG="Pushing the docs to $dir/ for branch: $CIRCLE_BRANCH, commit $CIRCLE_SHA1"

cd $HOME
if [ ! -d $DOC_REPO ];
then git clone --depth 1 --no-checkout "git@github.com:scikit-learn/"$DOC_REPO".git";
fi
cd $DOC_REPO

# check if it's a new branch

echo $dir > .git/info/sparse-checkout
if ! git show HEAD:$dir >/dev/null
then
	# directory does not exist. Need to make it so sparse checkout works
	mkdir $dir
	touch $dir/index.html
	git add $dir
fi
git checkout master
git reset --hard origin/master
if [ -d $dir ]
then
	git rm -rf $dir/ && rm -rf $dir/
fi
cp -R $GENERATED_DOC_DIR $dir
git config user.email "olivier.grisel+sklearn-ci@gmail.com"
git config user.name $USERNAME
git config push.default matching
git add -f $dir/
git commit -m "$MSG" $dir
git push
echo $MSG 
