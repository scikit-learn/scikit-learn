# This script is used by the CUDA CI workflow to comment on a PR
# as well as add a commit status to the commit it tests.
import os
from textwrap import dedent

import requests


def get_headers(token):
    """Get the headers for the GitHub API."""
    return {
        "Accept": "application/vnd.github+json",
        "Authorization": f"Bearer {token}",
        "X-GitHub-Api-Version": "2022-11-28",
    }


def find_first_bot_comment(repo, token, pr_number):
    """Get the first comment made by the bot."""
    # repo is in the form of "org/repo"
    # API doc: https://docs.github.com/en/rest/issues/comments?apiVersion=2022-11-28#list-issue-comments  # noqa
    response = requests.get(
        f"https://api.github.com/repos/{repo}/issues/{pr_number}/comments",
        headers=get_headers(token),
    )
    response.raise_for_status()
    all_comments = response.json()

    # Find all comments that match the CUDA CI bot, and return the first one.
    # There should always be only one such comment, or none, if the workflow
    # has not been run before
    comments = [
        comment
        for comment in all_comments
        if comment["user"]["login"] == "github-actions[bot]"
        and comment["body"].startswith("The CUDA CI bot")
    ]

    return comments[0] if comments else None


def create_bot_comment(message, repo, pr_number, token):
    """Create the initial comment"""
    response = requests.post(
        f"https://api.github.com/repos/{repo}/issues/{pr_number}/comments",
        headers=get_headers(token),
        json={"body": message},
    )

    response.raise_for_status()


def update_bot_comment(comment_id, message, repo, token):
    """Update the bot's comment"""
    response = requests.patch(
        f"https://api.github.com/repos/{repo}/issues/comments/{comment_id}",
        headers=get_headers(token),
        json={"body": message},
    )

    response.raise_for_status()


def add_commit_status(commit_id, workflow_status, details_url, token):
    # Translate from workflow step statuses to states that are available
    # for a commit status.
    if workflow_status == "success":
        state = "success"
    elif workflow_status == "skipped":
        # Unclear why the workflow status should ever be skipped, so
        # this probably means something went wrong
        state = "error"
    elif workflow_status == "cancelled":
        state = "pending"
    else:
        state = "failure"

    response = requests.post(
        f"https://api.github.com/repos/{repo}/statuses/{commit_id}",
        headers=get_headers(token),
        json={
            "state": state,
            "target_url": details_url,
            "description": "CUDA GPU tests",
            "context": "ci/cuda",
        },
    )

    response.raise_for_status()


if __name__ == "__main__":
    repo = os.environ["GITHUB_REPOSITORY"]
    token = os.environ["GITHUB_TOKEN"]
    pr_number = os.environ["PR_NUMBER"]
    sha = os.environ["BRANCH_SHA"]
    run_id = os.environ["RUN_ID"]
    workflow_status = os.environ["WORKFLOW_STATUS"]

    details_url = f"https://github.com/{repo}/actions/runs/{run_id}"

    comment = find_first_bot_comment(repo, token, pr_number)

    # end first line with \ to avoid the annoying empty line!
    message = f"""\
    The CUDA CI bot has run a workflow to check the changes of this PR ({sha}).

    Workflow status: {workflow_status}
    Details: {details_url}
    """
    message = dedent(message)

    if comment is None:
        create_bot_comment(message, repo, pr_number, token)
    else:
        update_bot_comment(comment["id"], message, repo, token)

    add_commit_status(sha, workflow_status, details_url, token)
