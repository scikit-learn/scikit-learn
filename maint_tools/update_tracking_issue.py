"""Creates or updates an issue if the CI fails. This is useful to keep track of
scheduled jobs that are failing repeatedly.

This script depends on:
- `defusedxml` for safer parsing for xml
- `PyGithub` for interacting with GitHub

The GitHub token only requires the `repo:public_repo` scope are described in
https://docs.github.com/en/developers/apps/building-oauth-apps/scopes-for-oauth-apps#available-scopes.
This scope allows the bot to create and edit its own issues. It is best to use a
github account that does **not** have commit access to the public repo.
"""

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import defusedxml.ElementTree as ET
from github import Github

parser = argparse.ArgumentParser(
    description="Create or update issue from JUnit test results from pytest"
)
parser.add_argument(
    "bot_github_token", help="Github token for creating or updating an issue"
)
parser.add_argument("ci_name", help="Name of CI run instance")
parser.add_argument("issue_repo", help="Repo to track issues")
parser.add_argument("link_to_ci_run", help="URL to link to")
parser.add_argument("--junit-file", help="JUnit file to determine if tests passed")
parser.add_argument(
    "--tests-passed",
    help=(
        "If --tests-passed is true, then the original issue is closed if the issue "
        "exists. If tests-passed is false, then the an issue is updated or created."
    ),
)
parser.add_argument(
    "--auto-close",
    help=(
        "If --auto-close is false, then issues will not auto close even if the tests"
        " pass."
    ),
    default="true",
)

args = parser.parse_args()

if args.junit_file is not None and args.tests_passed is not None:
    print("--junit-file and --test-passed can not be set together")
    sys.exit(1)

if args.junit_file is None and args.tests_passed is None:
    print("Either --junit-file or --test-passed must be passed in")
    sys.exit(1)

gh = Github(args.bot_github_token)
issue_repo = gh.get_repo(args.issue_repo)
dt_now = datetime.now(tz=timezone.utc)
date_str = dt_now.strftime("%b %d, %Y")
title_query = f"CI failed on {args.ci_name}"
title = f"⚠️ {title_query} (last failure: {date_str}) ⚠️"


def get_issue():
    login = gh.get_user().login
    issues = gh.search_issues(
        f"repo:{args.issue_repo} {title_query} in:title state:open author:{login}"
    )
    first_page = issues.get_page(0)
    # Return issue if it exist
    return first_page[0] if first_page else None


def create_or_update_issue(body=""):
    # Interact with GitHub API to create issue
    link = f"[{args.ci_name}]({args.link_to_ci_run})"
    issue = get_issue()

    max_body_length = 60_000
    original_body_length = len(body)
    # Avoid "body is too long (maximum is 65536 characters)" error from github REST API
    if original_body_length > max_body_length:
        body = (
            f"{body[:max_body_length]}\n...\n"
            f"Body was too long ({original_body_length} characters) and was shortened"
        )

    if issue is None:
        # Create new issue
        header = f"**CI failed on {link}** ({date_str})"
        issue = issue_repo.create_issue(title=title, body=f"{header}\n{body}")
        print(f"Created issue in {args.issue_repo}#{issue.number}")
        sys.exit()
    else:
        # Update existing issue
        header = f"**CI is still failing on {link}** ({date_str})"
        issue.edit(title=title, body=f"{header}\n{body}")
        print(f"Commented on issue: {args.issue_repo}#{issue.number}")
        sys.exit()


def close_issue_if_opened():
    print("Test has no failures!")
    issue = get_issue()
    if issue is not None:
        header_str = "## CI is no longer failing!"
        comment_str = (
            f"{header_str} ✅\n\n[Successful run]({args.link_to_ci_run}) on {date_str}"
        )

        print(f"Commented on issue #{issue.number}")
        # New comment if "## CI is no longer failing!" comment does not exist
        # If it does exist update the original comment which includes the new date
        for comment in issue.get_comments():
            if comment.body.startswith(header_str):
                comment.edit(body=comment_str)
                break
        else:  # no break
            issue.create_comment(body=comment_str)

        if args.auto_close.lower() == "true":
            print(f"Closing issue #{issue.number}")
            issue.edit(state="closed")
    sys.exit()


if args.tests_passed is not None:
    if args.tests_passed.lower() == "true":
        close_issue_if_opened()
    else:
        create_or_update_issue()

junit_path = Path(args.junit_file)
if not junit_path.exists():
    body = "Unable to find junit file. Please see link for details."
    create_or_update_issue(body)

# Find failures in junit file
tree = ET.parse(args.junit_file)
failure_cases = []

# Check if test collection failed
error = tree.find("./testsuite/testcase/error")
if error is not None:
    # Get information for test collection error
    failure_cases.append("Test Collection Failure")

for item in tree.iter("testcase"):
    failure = item.find("failure")
    if failure is None:
        continue

    failure_cases.append(item.attrib["name"])

if not failure_cases:
    close_issue_if_opened()

# Create content for issue
body_list = [f"- {case}" for case in failure_cases]
body = "\n".join(body_list)
create_or_update_issue(body)
