#!/bin/bash
# This script is meant to be called in the "deploy" step defined in 
# circle.yml. See https://circleci.com/docs/ for more details.
# The behavior of the script is controlled by environment variable defined
# in the circle.yml in the top level folder of the project.


if [ -z $CIRCLE_PROJECT_USERNAME ];
then USERNAME="sklearn-ci";
else USERNAME=$CIRCLE_PROJECT_USERNAME;
fi

DOC_REPO="scikit-learn.github.io"

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
git config core.sparseCheckout true
echo $dir > .git/info/sparse-checkout
git checkout $CIRCLE_BRANCH
git reset --hard origin/$CIRCLE_BRANCH
git rm -rf $dir/ && rm -rf $dir/
cp -R $HOME/scikit-learn/doc/_build/html/stable $dir
git config --global user.email "olivier.grisel+sklearn-ci@gmail.com"
git config --global user.name $USERNAME
git config --global push.default matching
git add -f $dir/
git commit -m "$MSG" $dir
git push
echo $MSG 

# List all available versions of the documentation
(echo "<html><body><h1>Available documentation for Scikit-learn</h1><ul>"
for d in $(git ls-tree --name-only master)
do
  # extract version number from Sphinx Javascript variable
  v=$(git cat-file blob master:$d/index.html 2>/dev/null | grep VERSION: | cut -d"'" -f2)
  if [ -n "$v" ]
  then
    echo "<li><a href=\"/$d\">Scikit-learn $v documentation</a>"
    if git cat-file -e master:$d/_downloads/scikit-learn-docs.pdf 2>/dev/null
    then
      echo "(<a href=\"/$d/_downloads/scikit-learn-docs.pdf\">PDF</a>)"
    fi
    echo "</li>"
  fi
done
echo "</ul></body></html>") > versions.html
git add versions.html
if git diff --cached --quiet
then
  MSG="Listing documentation versions in versions.html"
  git commit -m $MSG
  git push
fi
