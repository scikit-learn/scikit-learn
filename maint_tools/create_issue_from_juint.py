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

from pathlib import Path
import sys
import argparse

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
parser.add_argument("junit_file", help="JUnit file")
parser.add_argument(
    "max_failures",
    help="Maximum number of test failures to include in issue",
    default=5,
    type=int,
)

args = parser.parse_args()
gh = Github(args.bot_github_token)
issue_repo = gh.get_repo(args.issue_repo)
title = f"⚠️ CI failed on {args.ci_name} ⚠️"
max_failures = args.max_failures


def get_issue():
    login = gh.get_user().login
    issues = gh.search_issues(
        f"repo:{args.issue_repo} {title} in:title state:open author:{login}"
    )
    first_page = issues.get_page(0)
    # Return issue if it exist
    return first_page[0] if first_page else None


def create_or_update_issue(body):
    # Interact with GitHub API to create issue
    header = f"**CI Failed on [{args.ci_name}]({args.link_to_ci_run})**"
    body_text = f"{header}\n{body}"
    issue = get_issue()

    if issue is None:
        # Create new issue
        issue = issue_repo.create_issue(title=title, body=body_text)
        print(f"Created issue in {args.issue_repo}#{issue.number}")
        sys.exit()
    else:
        # Update existing issue
        issue.edit(title=title, body=body_text)
        print(f"Updated issue in {args.issue_repo}#{issue.number}")
        sys.exit()


junit_path = Path(args.junit_file)
if not junit_path.exists():
    body = "Unable to find junit file. Please see link for details."
    create_or_update_issue(body)
    sys.exit()

# Find failures in junit file
tree = ET.parse(args.junit_file)
failure_cases = []

# Check if test collection failed
error = tree.find("./testsuite/testcase/error")
if error is not None:
    # Get information for test collection error
    failure_cases.append({"title": "Test Collection Failure", "body": error.text})

for item in tree.iter("testcase"):
    failure = item.find("failure")
    if failure is None:
        continue

    failure_cases.append(
        {
            "title": item.attrib["name"],
            "body": failure.text,
        }
    )

if not failure_cases:
    print("Test has no failures!")
    issue = get_issue()
    if issue is not None:
        print(f"Closing issue #{issue.number}")
        new_body = (
            "## Closed issue because CI is no longer failing! ✅\n\n"
            f"[Successful run]({args.link_to_ci_run})\n\n"
            "## Previous failing issue\n\n"
            f"{issue.body}"
        )
        issue.edit(state="closed", body=new_body)
    sys.exit()

# Create content for issue
issue_summary = (
    "<details><summary>{title}</summary>\n\n```python\n{body}\n```\n</details>\n"
)
body_list = [issue_summary.format(**case) for case in failure_cases[:max_failures]]
body = "\n".join(body_list)
n_remaining_failures = len(failure_cases) - max_failures
if n_remaining_failures > 0:
    body += f"\n\nand [{n_remaining_failures} more failures]({args.link_to_ci_run})."
create_or_update_issue(body)
