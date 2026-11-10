"""Label contributor type using organization membership and merged PRs count.

Requires PyGithub and the GITHUB_TOKEN, READ_ORG_TOKEN, GITHUB_REPO and PR_NUMBER
environment variables. READ_ORG_TOKEN must be a PAT from a member of the
scikit-learn teams with read:org permission.

Pass --dry-run to report label changes without modifying the pull request.
"""

import argparse
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report label changes without modifying the pull request.",
    )
    args = parser.parse_args()
    repo_name = os.environ["GITHUB_REPO"]
    pr_number = int(os.environ["PR_NUMBER"])

    with (
        Github(auth=Auth.Token(os.environ["GITHUB_TOKEN"])) as gh,
        Github(auth=Auth.Token(os.environ["READ_ORG_TOKEN"])) as org_gh,
    ):
        repo = gh.get_repo(repo_name)
        pr = repo.get_pull(pr_number)
        author = pr.user.login
        org = org_gh.get_organization(repo_name.split("/", 1)[0])
        if org.has_in_members(pr.user):
            label = "contributor:team-member"
            print(f"{author}: member of {org.login} -> {label}")
        else:
            label = get_contributor_type(gh, repo_name, author)

        current_labels = {existing.name for existing in pr.get_labels()}

        if label not in current_labels:
            print(f"Adding {label}")
            if not args.dry_run:
                pr.add_to_labels(label)
        for candidate in CONTRIBUTOR_LABELS:
            if candidate != label and candidate in current_labels:
                print(f"Removing {candidate}")
                if not args.dry_run:
                    pr.remove_from_labels(candidate)


if __name__ == "__main__":
    main()
