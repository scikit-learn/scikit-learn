"""Label contributor type using the PR author's association and merged PRs count.

Requires PyGithub and the GITHUB_TOKEN, GITHUB_REPO and PR_NUMBER environment
variables.
"""

import os

from github import Auth, Github

CONTRIBUTOR_LABELS = (
    "contributor:first-time",
    "contributor:early",
    "contributor:seasoned",
    "contributor:team-member",
)


def get_contributor_type(gh: Github, repo_name: str, username: str) -> str:
    """Return a contributor label based on the user's merged PRs in repo_name."""
    # Read the raw response to access the incomplete_results flag.
    _, search = gh.requester.requestJsonAndCheck(
        "GET",
        "/search/issues",
        parameters={
            "q": f"repo:{repo_name} is:pr is:merged author:{username}",
            "per_page": 1,
        },
    )
    merged_count = search["total_count"]
    if search.get("incomplete_results"):
        print(
            "::warning::Merged-PR search returned incomplete results; "
            f"continuing with count {merged_count} for {username}."
        )
    if merged_count == 0:
        label = "contributor:first-time"
    elif merged_count < 10:
        label = "contributor:early"
    else:
        label = "contributor:seasoned"

    print(f"{username}: {merged_count} merged PRs in {repo_name} -> {label}")

    return label


def main():
    repo_name = os.environ["GITHUB_REPO"]
    pr_number = int(os.environ["PR_NUMBER"])

    with Github(auth=Auth.Token(os.environ["GITHUB_TOKEN"])) as gh:
        repo = gh.get_repo(repo_name)
        pr = repo.get_pull(pr_number)
        author = pr.user.login
        association = pr.author_association
        if association in {"MEMBER", "OWNER"}:
            label = "contributor:team-member"
            print(f"{author}: {association} -> {label}")
        else:
            label = get_contributor_type(gh, repo_name, author)

        current_labels = {existing.name for existing in pr.get_labels()}

        if label not in current_labels:
            print(f"Adding {label}")
            pr.add_to_labels(label)
        for candidate in CONTRIBUTOR_LABELS:
            if candidate != label and candidate in current_labels:
                print(f"Removing {candidate}")
                pr.remove_from_labels(candidate)


if __name__ == "__main__":
    main()
