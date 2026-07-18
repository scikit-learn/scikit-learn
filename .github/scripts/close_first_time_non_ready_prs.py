"""Close PRs from first-time contributors that link to a not-ready issue.

Called from .github/workflows/close-first-time-non-ready-prs.yml.
"""

import json
import os
import urllib.request

from github import Github

GITHUB_REPO = "scikit-learn/scikit-learn"
GITHUB_TOKEN = os.getenv("GITHUB_TOKEN")
PR_NUMBER = int(os.getenv("PR_NUMBER"))

CLOSING_ISSUE_QUERY = """
query($owner: String!, $name: String!, $pr_number: Int!) {
  repository(owner: $owner, name: $name) {
    pullRequest(number: $pr_number) {
      closingIssuesReferences(first: 1) {
        nodes {
          number
          labels(first: 20) {
            nodes {
              name
            }
          }
        }
      }
    }
  }
}
"""


def get_linked_issue(gh_repo, pr_number):
    owner, name = gh_repo.split("/")
    request = urllib.request.Request(
        "https://api.github.com/graphql",
        method="POST",
        headers={
            "Authorization": f"Bearer {GITHUB_TOKEN}",
            "Content-Type": "application/json",
        },
        data=json.dumps(
            {
                "query": CLOSING_ISSUE_QUERY,
                "variables": {"owner": owner, "name": name, "pr_number": pr_number},
            }
        ).encode(),
    )
    with urllib.request.urlopen(request) as response:
        payload = json.load(response)
    nodes = payload["data"]["repository"]["pullRequest"]["closingIssuesReferences"][
        "nodes"
    ]
    return nodes[0] if nodes else None


def is_not_ready(issue):
    return any(
        label["name"].startswith("Needs") or label["name"] == "RFC"
        for label in issue["labels"]["nodes"]
    )


linked_issue = get_linked_issue(GITHUB_REPO, PR_NUMBER)

if linked_issue and is_not_ready(linked_issue):
    gh = Github(GITHUB_TOKEN)
    repo = gh.get_repo(GITHUB_REPO)
    pr = repo.get_pull(PR_NUMBER)

    MESSAGE = (
        "Thank you for your interest in contributing to scikit-learn.\n\n"
        "The linked issue is still under discussion, and the maintainers have not "
        "yet reached consensus on how it should be resolved. Before opening a pull "
        "request, please review the \"Issues tagged 'Needs Triage'\" section of the "
        "Contributing Guide:\n"
        "https://scikit-learn.org/dev/developers/contributing.html"
        "#issues-tagged-needs-triage\n\n"
        "For now, we are closing this pull request.\n\n"
        "If you believe your proposed change addresses the linked issue, please "
        "explain your proposal in that issue, and allow time for discussion with "
        "the maintainers so that consensus can be reached before implementation "
        "begins."
    )

    print(f"Closing PR #{PR_NUMBER} with comment")
    pr.create_issue_comment(MESSAGE)
    pr.edit(state="closed")
