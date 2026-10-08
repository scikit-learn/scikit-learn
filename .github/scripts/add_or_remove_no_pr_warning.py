# Used in Github Action .github/workflows/not_ready_for_pr_warning.yml
import argparse
import os

from github import Auth, Github

parser = argparse.ArgumentParser(
    description="Add or remove no-contribution warning from issue body"
)
parser.add_argument(
    "--mode", choices=["add", "remove"], help="Whether to add or remove warning"
)
args = parser.parse_args()

# env variables are defined in .github/workflows/not_ready_for_pr_warning.yml
g = Github(auth=Auth.Token(os.environ["GITHUB_TOKEN"]))
repo = g.get_repo(os.environ["GITHUB_REPO"])
issue = repo.get_issue(number=int(os.environ["ISSUE_NUMBER"]))

body_text = str(issue.body) if issue.body else ""

message = (
    "> [!WARNING]\n"
    "> This issue is not yet ready for a PR. If you are interested in contributing to "
    "scikit-learn, please have a look at our [contributing guidelines]"
    "(https://scikit-learn.org/dev/developers/contributing.html), and in particular "
    "the sections for [new contributors]"
    "(https://scikit-learn.org/dev/developers/contributing.html#new-contributors) and "
    'the ["Needs triage"](https://scikit-learn.org/dev/developers/contributing.html#'
    "issues-tagged-needs-triage) label."
)

if args.mode == "add":
    if not body_text.startswith(message):
        new_body = f"{message}\n\n{body_text}"
        issue.edit(body=new_body)
        print(f"Added warning to issue: {os.environ['GITHUB_REPO']}#{issue.number}")

else:
    still_needs_something = any(
        label.name.startswith("Needs") or label.name == "RFC" for label in issue.labels
    )
    if not still_needs_something:
        if body_text.startswith(message):
            new_body = body_text.removeprefix(f"{message}\n\n")
            issue.edit(body=new_body)
            print(
                "Removed warning from issue: "
                f"{os.environ['GITHUB_REPO']}#{issue.number}"
            )
