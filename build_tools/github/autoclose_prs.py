"""Close PRs labeled with 'autoclose' more than 14 days ago.

Called from .github/workflows/autoclose-schedule.yml."""

import os
from datetime import datetime, timedelta, timezone
from pprint import pprint

from github import Auth, Github


def get_labeled_last_time(pr, label):
    labeled_time = datetime.max
    for event in pr.get_events():
        if event.event == "labeled" and event.label.name == label:
            labeled_time = event.created_at

    return labeled_time


dry_run = False
cutoff_days = 14

gh_repo = "scikit-learn/scikit-learn"
github_token = os.getenv("GITHUB_TOKEN")

auth = Auth.Token(github_token)
gh = Github(auth=auth)
repo = gh.get_repo(gh_repo)


now = datetime.now(timezone.utc)
label = "autoclose"
prs = [
    each for each in repo.get_issues(labels=[label]) if each.pull_request is not None
]
prs_info = [f"{pr.title}: {pr.html_url}" for pr in prs]
print(f"Found {len(prs)} opened PRs with label {label}")
pprint(prs_info)

prs = [
    pr
    for pr in prs
    if (now - get_labeled_last_time(pr, label)) > timedelta(days=cutoff_days)
]
prs_info = [f"{pr.title} {pr.html_url}" for pr in prs]
print(f"Found {len(prs)} PRs to autoclose")
pprint(prs_info)

message = (
    "Thank you for your interest in contributing to scikit-learn, but we cannot "
    "accept your contribution as this pull request does not meet our development "
    "standards.\n\n"
    "Following our autoclose policy, we are closing this PR after allowing two "
    "weeks time for improvements.\n\n"
    "Thank you for your understanding. If you think your PR has been closed "
    "by mistake, please comment below."
)

for pr in prs:
    print(f"Closing PR #{pr.number} with comment")
    if not dry_run:
        pr.create_comment(message)
        pr.edit(state="closed")
