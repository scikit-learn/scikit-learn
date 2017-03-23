# Get current commit hash
export COMMIT_HASH=$(git log -n 1 --pretty=format:"%H")
# Get the git root directory
export sklearn_root=$(git rev-parse --show-toplevel)

affected_doc_paths() {
	files=$(git diff --name-only origin/master...$COMMIT_HASH)
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
echo $affected
(
echo '<html><body><ul>'
echo "$affected" | sed 's|.*|<li><a href="&">&</a></li>|'
echo '</ul></body></html>'
) > $sklearn_root'/doc/_build/html/stable/_changed.html'
