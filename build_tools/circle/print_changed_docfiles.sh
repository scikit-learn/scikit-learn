#!/usr/bin/env bash
set -e

# Create a _changed.html file inside the build directory 
# in "doc/_build/html/stable".
# This file contains the list of the files that may have changed 
# due to the modification of code.
#  
# The documentation needs to be built before this script is executed.
#
# The script can be either launched directly from the terminal:
# $> source build_tools/circle/print_changed_docfiles.sh
# Or from the doc folder with a make command
# $> make changed-from-master
#
# The script accepts a parameter, the git reference base, if no value 
# is provided, the ref base is set to origin/master.

# Get current commit hash
COMMIT_HASH=$(git log -n 1 --pretty=format:"%H")
# Get the git root directory
sklearn_root=$(git rev-parse --show-toplevel)

if [ ! -d "$sklearn_root/doc/_build/html/stable" ]
then
    echo "The folder \"doc/_build/html/stable\" does not exist."
    echo "You need to build the documentation first."
    exit
fi

# Parameter 1: the git ref base (Default: origin/master)
if [ ! -n "$1" ]
then
    ref_base=origin/master
else
    ref_base=$1
fi

affected_doc_paths() {
	files=$(git diff --name-only $ref_base...$COMMIT_HASH)
	echo "$files" | grep ^doc/.*\.rst | sed 's/^doc\/\(.*\)\.rst$/\1.html/'
	echo "$files" | grep ^examples/.*.py | sed 's/^\(.*\)\.py$/auto_\1.html/'
	sklearn_files=$(echo "$files" | grep '^sklearn/')
	if [ -n "$sklearn_files" ]
	then
		grep -hlR -f<(echo "$sklearn_files" | sed 's/^/scikit-learn\/blob\/[a-z0-9]*\//') $sklearn_root/doc/_build/html/stable/modules/generated | awk -F"stable/" '{print $2}'
	fi
}

echo "The following documentation files may have been changed by commit $COMMIT_HASH:"
affected=$(affected_doc_paths)
if [ -n "$affected" ]
then
    echo "Files modificated found and written in _build/html/stable/_changed.html"
    (
    echo '<html><body><ul>'
    echo "$affected" | sed 's|.*|<li><a href="&">&</a></li>|'
    echo '</ul></body></html>'
    ) > $sklearn_root'/doc/_build/html/stable/_changed.html'
else
    echo 'No file modificated found in the documentation.'
    (
    echo '<html><body>'
    echo 'No file modificated found in the documentation.'
    echo '</body></html>'
    ) > $sklearn_root'/doc/_build/html/stable/_changed.html'
fi
