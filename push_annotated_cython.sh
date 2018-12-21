#!/bin/sh

set -e  # exit if any command fails


BRANCH=gbm
SOURCE_DIR=/home/nico/dev/sklearn/sklearn/ensemble/gbm
TARGET_DIR=/home/nico/dev/cython_annotations

ORIGINAL_DIR=`pwd`


git co $BRANCH

# Commits in the branch (provided it branched off master)
COMMITS=`git log master.. --pretty=format:"%h"`

annotate_and_copy_files() {
  # For a give commit, annotate all pyx file in SOURCE_DIR and copy the html
  # files in TARGET_DIR/COMMIT_HASH/

  git co $1  # checkout commit
  for pyx_file in `ls $SOURCE_DIR/*.pyx`
  do
    echo 'annotating' $1 $pyx_file
    cython -a $pyx_file
  done

  for html_file in `ls $SOURCE_DIR/*.html`
  do
    mkdir -p $TARGET_DIR/$1
    cp $html_file $TARGET_DIR/$1
    html_file_name=$(basename -- "$html_file")  # without path
    echo Copied $html_file_name to $TARGET_DIR/$1
  done
}

for commit in $COMMITS
do
  annotate_and_copy_files $commit
done


# Get into target dir, commit html files and push them.
cd $TARGET_DIR
git co gh-pages
echo Generating index.html
python lol.py  # generates index.html with links to each file
echo Committing and pushing files
git add .
git ci -am "Added some annotated cython files"
git push

cd $ORIGINAL_DIR  # go back where we were
git co $BRANCH  # Probably useless since with checked out the last commit
