"""Labels PRs based on title. Must be run in a github action with the
pull_request_target event."""
from ghapi.all import context_github
from ghapi.all import GhApi
from ghapi.all import user_repo
from ghapi.all import github_token
import re

owner, repo = user_repo()
pull_request = context_github.event.pull_request
title = pull_request.title

regex_to_labels = [
    (r"\bDOC\b", "Documentation"),
    (r"\bCI\b", "Build / CI")
]

labels_to_add = [
    label for regex, label in regex_to_labels
    if re.search(regex, title)
]

if labels_to_add:
    api = GhApi(owner=owner, repo=repo, token=github_token())
    api.issues.add_labels(pull_request.number, labels=labels_to_add)
