"""First-time Contributor Workflow

Called from .github/workflows/welcome-first-time-contributor.yml.
"""

import json
import os
import urllib.request

from github import Github

GITHUB_REPO = os.getenv("GITHUB_REPO")
GITHUB_TOKEN = os.getenv("GITHUB_TOKEN")
PR_NUMBER = int(os.getenv("PR_NUMBER"))


def get_linked_issue(gh_repo, pr_number):
    """Get the linked issue of the PR."""
    owner, name = gh_repo.split("/")

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
    """Determine if an issue is still not ready for a PR."""
    return any(
        label["name"].startswith("Needs") or label["name"] == "RFC"
        for label in issue["labels"]["nodes"]
    )


gh = Github(GITHUB_TOKEN)
repo = gh.get_repo(GITHUB_REPO)
pr = repo.get_pull(PR_NUMBER)

linked_issue = get_linked_issue(GITHUB_REPO, PR_NUMBER)

if linked_issue and is_not_ready(linked_issue):
    # Close the PR if the linked issue is not ready
    MESSAGE = (
        "Thank you for your interest in contributing to scikit-learn.\n\n"
        "The linked issue is still under discussion, and the maintainers have not "
        "yet reached consensus on how it should be resolved. Before opening a pull "
        "request, please review the \"Issues tagged 'Needs Triage'\" section of the "
        "Contributing Guide:\n"
        "https://scikit-learn.org/stable/developers/contributing.html"
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
    pr.add_to_labels("autoclose")

if pr.state == "open":
    # Post welcome comment
    MESSAGE = (
        "Thank you for opening your first pull request to scikit-learn! 🎉"
        "\n\n"
        "To help get your contribution reviewed, please make sure that:"
        "\n\n"
        "* You have filled out the "
        "[pull request template]"
        "(https://github.com/scikit-learn/scikit-learn/blob/main/.github/PULL_REQUEST_TEMPLATE.md)."
        "\n\n"
        "* The pull request addresses an existing issue that is ready for contribution "
        "(e.g. not tagged as 'Needs Triage', 'Needs Decision', ...). "
        "If you are proposing a new feature, please open an issue to discuss it first."
        "\n\n"
        "* There are no other open pull requests already targeting the same issue."
        "\n\n"
        "* You have followed the "
        "[pull request checklist]"
        "(https://scikit-learn.org/stable/developers/contributing.html#pull-request-checklist)."
        " In particular, linting and tests should pass."
    )

    pr.create_issue_comment(MESSAGE)
