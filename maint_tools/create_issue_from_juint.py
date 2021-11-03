"""Creates or updates an issue if the CI fails. This is useful to keep track of
scheduled jobs that an fail.

This scirpt depends on:
- defusedxml for safer parsing for xml
- PyGithub for interacting with GitHub
"""

from pathlib import Path
import sys
import argparse

import defusedxml.ElementTree as ET
from github import Github

# Labels to place on issue
CI_LABEL_NAMES = ["Build / CI"]

parser = argparse.ArgumentParser(
    description="Create or update issue from JUnit test results from pytest"
)
parser.add_argument(
    "bot_github_token", help="Github token for creating or updating an issue"
)
parser.add_argument("ci_name", help="Name of CI run instance")
parser.add_argument("issue_repo", help="Repo to track issues")
parser.add_argument("link_to_run", help="URL to link to")
parser.add_argument("junit_file", help="JUnit file")

args = parser.parse_args()


def create_or_update_issue(body):
    # Interact with GitHub API to create issue
    title = f"⚠️ CI failed on {args.ci_name} ⚠️"
    header = f"**CI Failed on [{args.ci_name}]({args.link_to_run})**"
    body_text = f"{header}\n{body}"

    gh = Github(args.bot_github_token)
    issue_repo = gh.get_repo(args.issue_repo)
    login = gh.get_user().login
    issues = gh.search_issues(
        f"repo:{args.issue_repo} {title} in:title state:open author:{login}"
    )

    first_page = issues.get_page(0)
    if not first_page:
        # Create new issue
        labels = [issue_repo.get_label(label) for label in CI_LABEL_NAMES]
        issue = issue_repo.create_issue(title, body=body_text, labels=labels)
        print(f"Created issue in {args.issue_repo}#{issue.number}")
        sys.exit(0)
    else:
        # Update existing issue
        issue = first_page[0]
        issue.edit(title, body=body_text)
        print(f"Updated issue in {args.issue_repo}#{issue.number}")
        sys.exit(0)


junit_path = Path(args.junit_file)
if not junit_path.exists():
    body = "Unable to find junit file. Please see link to details."
    create_or_update_issue(body)
    sys.exit(1)

# Find failures in junit file
tree = ET.parse(args.junit_file)
failure_cases = []

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
    # TODO: Do we want to close the issue here?
    print("Test has no failures!")
    sys.exit(0)

# Create content for issue
issue_summary = (
    "<details><summary>{title}</summary>\n\n```python\n{body}\n```\n</details>\n"
)
body_list = [issue_summary.format(**case) for case in failure_cases]
body = "\n".join(body_list)
create_or_update_issue(body)
