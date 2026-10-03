"""Welcome first-time contributors.

Called from .github/workflows/first-time-contributor.yml.
"""

import os

from github import Auth, Github

GITHUB_REPO = os.getenv("GITHUB_REPO")
GITHUB_TOKEN = os.getenv("GITHUB_TOKEN")
PR_NUMBER = int(os.getenv("PR_NUMBER"))

gh = Github(auth=Auth.Token(GITHUB_TOKEN))
repo = gh.get_repo(GITHUB_REPO)
pr = repo.get_pull(PR_NUMBER)

if pr.state == "open" and pr.author_association in {
    "FIRST_TIME_CONTRIBUTOR",
    "FIRST_TIMER",
    "NONE",
}:
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

    print(f"Posting welcome comment in #{PR_NUMBER}")
    pr.create_issue_comment(MESSAGE)
