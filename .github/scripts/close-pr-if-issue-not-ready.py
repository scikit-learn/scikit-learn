"""Close PRs linked to unready issues for early contributors.

Early contributors have fewer than 10 merged PRs in this repository.

Called from .github/workflows/close-pr-if-issue-not-ready.yml.
"""

import os

from github import Auth, Github

GITHUB_REPO = os.getenv("GITHUB_REPO")
GITHUB_TOKEN = os.getenv("GITHUB_TOKEN")
PR_NUMBER = int(os.getenv("PR_NUMBER"))


def get_linked_issue(gh, gh_repo, pr_number):
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
    _, payload = gh.requester.graphql_query(
        CLOSING_ISSUE_QUERY,
        {"owner": owner, "name": name, "pr_number": pr_number},
    )
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


gh = Github(auth=Auth.Token(GITHUB_TOKEN))
repo = gh.get_repo(GITHUB_REPO)
pr = repo.get_pull(PR_NUMBER)

linked_issue = get_linked_issue(gh, GITHUB_REPO, PR_NUMBER)

if linked_issue and is_not_ready(linked_issue):
    merged_pr_count = gh.search_issues(
        f"repo:{GITHUB_REPO} is:pr is:merged author:{pr.user.login}"
    ).totalCount
else:
    merged_pr_count = None

if merged_pr_count is not None and merged_pr_count < 10:
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
        "* If you believe your proposed change addresses the linked issue, please "
        "explain your proposal in that issue, and allow time for discussion with "
        "the maintainers so that consensus can be reached before implementation "
        "begins.\n\n"
        "* If you believe this is a mistake, please leave a comment on "
        "the linked issue as well."
    )

    print(f"Closing PR #{PR_NUMBER} with comment")
    pr.create_issue_comment(MESSAGE)
    pr.edit(state="closed")
    pr.add_to_labels("linked issue not ready")
