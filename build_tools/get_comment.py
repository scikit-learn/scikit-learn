# This script is used to generate a comment for a PR when linting issues are
# detected. It is used by the `Comment on failed linting` GitHub Action.
# This script fails if there are not comments to be posted.

import os
import re

import requests


def get_versions(versions_file):
    """Get the versions of the packages used in the linter job.

    Parameters
    ----------
    versions_file : str
        The path to the file that contains the versions of the packages.

    Returns
    -------
    versions : dict
        A dictionary with the versions of the packages.
    """
    with open(versions_file, "r") as f:
        return dict(line.strip().split("=") for line in f)


def get_step_message(log, start, end, title, message, details):
    """Get the message for a specific test.

    Parameters
    ----------
    log : str
        The log of the linting job.

    start : str
        The string that marks the start of the test.

    end : str
        The string that marks the end of the test.

    title : str
        The title for this section.

    message : str
        The message to be added at the beginning of the section.

    details : bool
        Whether to add the details of each step.

    Returns
    -------
    message : str
        The message to be added to the comment.
    """
    if end not in log:
        return ""
    res = (
        f"-----------------------------------------------\n### {title}\n\n{message}\n\n"
    )
    if details:
        res += (
            "<details>\n\n```\n"
            + log[log.find(start) + len(start) + 1 : log.find(end) - 1]
            + "\n```\n\n</details>\n\n"
        )
    return res


def get_message(log_file, repo, pr_number, sha, run_id, details, versions):
    with open(log_file, "r") as f:
        log = f.read()

    sub_text = (
        "\n\n<sub> _Generated for commit:"
        f" [{sha[:7]}](https://github.com/{repo}/pull/{pr_number}/commits/{sha}). "
        "Link to the linter CI: [here]"
        f"(https://github.com/{repo}/actions/runs/{run_id})_ </sub>"
    )

    if "### Linting completed ###" not in log:
        return (
            "## ❌ Linting issues\n\n"
            "There was an issue running the linter job. Please update with "
            "`upstream/main` ([link]("
            "https://scikit-learn.org/dev/developers/contributing.html"
            "#how-to-contribute)) and push the changes. If you already have done "
            "that, please send an empty commit with `git commit --allow-empty` "
            "and push the changes to trigger the CI.\n\n" + sub_text
        )

    message = ""

    # ruff check
    message += get_step_message(
        log,
        start="### Running the ruff linter ###",
        end="Problems detected by ruff check",
        title="`ruff check`",
        message=(
            "`ruff` detected issues. Please run "
            "`ruff check --fix --output-format=full` locally, fix the remaining "
            "issues, and push the changes. Here you can see the detected issues. Note "
            f"that the installed `ruff` version is `ruff={versions['ruff']}`."
        ),
        details=details,
    )

    # ruff format
    message += get_step_message(
        log,
        start="### Running the ruff formatter ###",
        end="Problems detected by ruff format",
        title="`ruff format`",
        message=(
            "`ruff` detected issues. Please run `ruff format` locally and push "
            "the changes. Here you can see the detected issues. Note that the "
            f"installed `ruff` version is `ruff={versions['ruff']}`."
        ),
        details=details,
    )

    # mypy
    message += get_step_message(
        log,
        start="### Running mypy ###",
        end="Problems detected by mypy",
        title="`mypy`",
        message=(
            "`mypy` detected issues. Please fix them locally and push the changes. "
            "Here you can see the detected issues. Note that the installed `mypy` "
            f"version is `mypy={versions['mypy']}`."
        ),
        details=details,
    )

    # cython-lint
    message += get_step_message(
        log,
        start="### Running cython-lint ###",
        end="Problems detected by cython-lint",
        title="`cython-lint`",
        message=(
            "`cython-lint` detected issues. Please fix them locally and push "
            "the changes. Here you can see the detected issues. Note that the "
            "installed `cython-lint` version is "
            f"`cython-lint={versions['cython-lint']}`."
        ),
        details=details,
    )

    # deprecation order
    message += get_step_message(
        log,
        start="### Checking for bad deprecation order ###",
        end="Problems detected by deprecation order check",
        title="Deprecation Order",
        message=(
            "Deprecation order check detected issues. Please fix them locally and "
            "push the changes. Here you can see the detected issues."
        ),
        details=details,
    )

    # doctest directives
    message += get_step_message(
        log,
        start="### Checking for default doctest directives ###",
        end="Problems detected by doctest directive check",
        title="Doctest Directives",
        message=(
            "doctest directive check detected issues. Please fix them locally and "
            "push the changes. Here you can see the detected issues."
        ),
        details=details,
    )

    # joblib imports
    message += get_step_message(
        log,
        start="### Checking for joblib imports ###",
        end="Problems detected by joblib import check",
        title="Joblib Imports",
        message=(
            "`joblib` import check detected issues. Please fix them locally and "
            "push the changes. Here you can see the detected issues."
        ),
        details=details,
    )

    if not message:
        # no issues detected, so this script "fails"
        return (
            "## ✔️ Linting Passed\n"
            "All linting checks passed. Your pull request is in excellent shape! ☀️"
            + sub_text
        )

    if not details:
        # This happens if posting the log fails, which happens if the log is too
        # long. Typically, this happens if the PR branch hasn't been updated
        # since we've introduced import sorting.
        branch_not_updated = (
            "_Merging with `upstream/main` might fix / improve the issues if you "
            "haven't done that since 21.06.2023._\n\n"
        )
    else:
        branch_not_updated = ""

    message = (
        "## ❌ Linting issues\n\n"
        + branch_not_updated
        + "This PR is introducing linting issues. Here's a summary of the issues. "
        + "Note that you can avoid having linting issues by enabling `pre-commit` "
        + "hooks. Instructions to enable them can be found [here]("
        + "https://scikit-learn.org/dev/developers/contributing.html#how-to-contribute)"
        + ".\n\n"
        + "You can see the details of the linting issues under the `lint` job [here]"
        + f"(https://github.com/{repo}/actions/runs/{run_id})\n\n"
        + message
        + sub_text
    )

    return message


def get_headers(token):
    """Get the headers for the GitHub API."""
    return {
        "Accept": "application/vnd.github+json",
        "Authorization": f"Bearer {token}",
        "X-GitHub-Api-Version": "2022-11-28",
    }


def find_lint_bot_comments(repo, token, pr_number):
    """Get the comment from the linting bot."""
    # repo is in the form of "org/repo"
    # API doc: https://docs.github.com/en/rest/issues/comments?apiVersion=2022-11-28#list-issue-comments
    response = requests.get(
        f"https://api.github.com/repos/{repo}/issues/{pr_number}/comments",
        headers=get_headers(token),
    )
    response.raise_for_status()
    all_comments = response.json()

    failed_comment = "❌ Linting issues"
    success_comment = "✔️ Linting Passed"

    # Find all comments that match the linting bot, and return the first one.
    # There should always be only one such comment, or none, if the PR is
    # just created.
    comments = [
        comment
        for comment in all_comments
        if comment["user"]["login"] == "github-actions[bot]"
        and (failed_comment in comment["body"] or success_comment in comment["body"])
    ]

    if len(all_comments) > 25 and not comments:
        # By default the API returns the first 30 comments. If we can't find the
        # comment created by the bot in those, then we raise and we skip creating
        # a comment in the first place.
        raise RuntimeError("Comment not found in the first 30 comments.")

    return comments[0] if comments else None


def create_or_update_comment(comment, message, repo, pr_number, token):
    """Create a new comment or update existing one."""
    # repo is in the form of "org/repo"
    if comment is not None:
        print("updating existing comment")
        # API doc: https://docs.github.com/en/rest/issues/comments?apiVersion=2022-11-28#update-an-issue-comment
        response = requests.patch(
            f"https://api.github.com/repos/{repo}/issues/comments/{comment['id']}",
            headers=get_headers(token),
            json={"body": message},
        )
    else:
        print("creating new comment")
        # API doc: https://docs.github.com/en/rest/issues/comments?apiVersion=2022-11-28#create-an-issue-comment
        response = requests.post(
            f"https://api.github.com/repos/{repo}/issues/{pr_number}/comments",
            headers=get_headers(token),
            json={"body": message},
        )

    response.raise_for_status()


if __name__ == "__main__":
    repo = os.environ["GITHUB_REPOSITORY"]
    token = os.environ["GITHUB_TOKEN"]
    pr_number = os.environ["PR_NUMBER"]
    sha = os.environ["BRANCH_SHA"]
    log_file = os.environ["LOG_FILE"]
    run_id = os.environ["RUN_ID"]
    versions_file = os.environ["VERSIONS_FILE"]

    versions = get_versions(versions_file)

    if not repo or not token or not pr_number or not log_file or not run_id:
        raise ValueError(
            "One of the following environment variables is not set: "
            "GITHUB_REPOSITORY, GITHUB_TOKEN, PR_NUMBER, LOG_FILE, RUN_ID"
        )

    if not re.match(r"\d+$", pr_number):
        raise ValueError(f"PR_NUMBER should be a number, got {pr_number!r} instead")

    try:
        comment = find_lint_bot_comments(repo, token, pr_number)
    except RuntimeError:
        print("Comment not found in the first 30 comments. Skipping!")
        exit(0)

    try:
        message = get_message(
            log_file,
            repo=repo,
            pr_number=pr_number,
            sha=sha,
            run_id=run_id,
            details=True,
            versions=versions,
        )
        create_or_update_comment(
            comment=comment,
            message=message,
            repo=repo,
            pr_number=pr_number,
            token=token,
        )
        print(message)
    except requests.HTTPError:
        # The above fails if the message is too long. In that case, we
        # try again without the details.
        message = get_message(
            log_file,
            repo=repo,
            pr_number=pr_number,
            sha=sha,
            run_id=run_id,
            details=False,
            versions=versions,
        )
        create_or_update_comment(
            comment=comment,
            message=message,
            repo=repo,
            pr_number=pr_number,
            token=token,
        )
        print(message)
