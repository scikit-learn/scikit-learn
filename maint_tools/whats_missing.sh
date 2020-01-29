#!/bin/bash
# This script helps identify pull requests that were merged without a what's
# new entry, where one would be appropriate.

if [ $# -ne 2 ]
then
	echo "Usage: GITHUB_TOKEN=... $0 <prev_release_ref> <whats_new_version>" >&2
	exit 1
fi
from_branch=$1
to_file=$2

logged_prs() {
	git log --oneline $from_branch..master sklearn/ |
		grep -wv -e CLN -e TST -e CI -e DOC -e doc -e MNT -e MAINT -e BLD -e COSMIT -e EXA -e examples -e example -e minor -e STY -e Style -e docstring |
		grep -o '(#[0-9][0-9]\+)$' |
		grep -o '[0-9]\+'
}

mentioned_issues() {
	cat doc/whats_new/v$to_file.rst |
			grep -o 'issue:`[0-9]\+`' |
			grep -o '[0-9]\+'
}

get_closed_issues() {
	pr=$1
	url=https://api.github.com/repos/scikit-learn/scikit-learn/pulls/$pr
	python - $url <<EOF
import json
import sys
import re
import os
from urllib import request

req = request.Request(sys.argv[1], headers={"Authorization": "token %s" % os.environ['GITHUB_TOKEN']})
body = json.loads(request.urlopen(req).read().decode('utf8'))['body']
body = re.sub('<!--.*?-->', '', body, flags=re.DOTALL)
matches = re.findall(r'(?i)\\b(?:fix|fixes|resolve|resolves|close|closes) +(?:https?://github.com/scikit-learn/scikit-learn/(?:pull|issues)/|#)?([0-9]+)',
                          body)
print(' '.join(matches))
EOF
}

pr_numbers=$(diff <(logged_prs | sort) <(mentioned_issues | sort) |
	grep '<' |
	cut -c3- |
	grep -v -w -Ff <(git log --oneline $from_branch | grep -o '(#[0-9][0-9]\+)$' | grep -o '[0-9]\+') )  # drop things already released

filtered_pr_numbers=$(
	for pr in $pr_numbers
	do
		echo $pr $(get_closed_issues $pr)
	done |
		grep -v -wFf <(mentioned_issues) |
		cut -d' ' -f1
)

echo $filtered_pr_numbers |
	sed 's/[^ ]*/--grep (#&)/g' |
	xargs git log
