"""Close PRs labeled with 'autoclose' more than 14 days ago.

Called from .github/workflows/autoclose-schedule.yml."""

import os
from datetime import datetime, timezone

from github import Github

CUTOFF_DAYS = 14


def get_labeled_last_time(pr, label):
    labeled_time = None
    for event in pr.get_events():
        if event.event == "labeled" and event.label.name == label:
            labeled_time = event.created_at

    return labeled_time


gh_repo = "scikit-learn/scikit-learn"
github_token = os.getenv("GITHUB_TOKEN")

gh = Github(github_token)
repo = gh.get_repo(gh_repo)


now = datetime.now(timezone.utc)
label = "autoclose"
prs = [
    each
    for each in repo.get_issues(labels=[label])
    if each.pull_request is not None
    and (now - get_labeled_last_time(each, label)).days > CUTOFF_DAYS
]
pr_numbers = [pr.number for pr in prs]
print(f"Found {len(prs)} PRs to autoclose: {pr_numbers}")

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
    pr.create_comment(message)
    pr.edit(state="closed")
