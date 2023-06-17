# This script is used to generate a comment for a PR when linting issues are
# detected. It is used by the `Comment on failed linting` GitHub Action.
# This script fails if there are not comments to be posted.

import os
import requests


def get_step_message(log, start, end, title, message):
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

    Returns
    -------
    message : str
        The message to be added to the comment.
    """
    if end not in log:
        return ""
    return (
        "-----------------------------------------------\n"
        + f"### {title}\n\n"
        + message
        + "\n\n<details>\n\n```\n"
        + log[log.find(start) + len(start) + 1 : log.find(end) - 1]
        + "\n```\n\n</details>\n\n"
    )


def get_message():
    with open("linting_output.txt", "r") as f:
        log = f.read()

    message = ""

    # black
    message += get_step_message(
        log,
        start="### Running black ###",
        end="Problems detected by black",
        title="`black`",
        message=(
            "`black` detected issues. Please run `black .` locally and push "
            "the changes. Here you can see the detected issues. Note that "
            "running black might also fix some of the issues which might be "
            "detected by `flake8`."
        ),
    )

    # flake8
    message += get_step_message(
        log,
        start="### Running flake8 ###",
        end="Problems detected by flake8",
        title="`flake8`",
        message=(
            "`flake8` detected issues. Please fix them locally and push the changes. "
            "Here you can see the detected issues."
        ),
    )

    # mypy
    message += get_step_message(
        log,
        start="### Running mypy ###",
        end="Problems detected by mypy",
        title="`mypy`",
        message=(
            "`mypy` detected issues. Please fix them locally and push the changes. "
            "Here you can see the detected issues."
        ),
    )

    # cython-lint
    message += get_step_message(
        log,
        start="### Running cython-lint ###",
        end="Problems detected by cython-lint",
        title="`cython-lint`",
        message=(
            "`cython-lint` detected issues. Please fix them locally and push "
            "the changes. Here you can see the detected issues."
        ),
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
    )

    if not len(message):
        # no issues detected, so this script "fails"
        return (
            "## Linting Passed\n"
            "All linting checks passed. Your pull request is in excellent shape!"
        )

    message = (
        "## Linting issues\n\n"
        "This PR is introducing linting issues. Here's a summary of the issues. "
        "Note that you can avoid having linting issues by enabling `pre-commit` "
        "hooks. Instructions to enable them can be found [here]("
        "https://scikit-learn.org/dev/developers/contributing.html#how-to-contribute)."
        "\n\n"
        + message
    )

    return message


def get_headers(token):
    """Get the headers for the GitHub API."""
    return {
        "Accept": "application/vnd.github+json",
        "Authorization": f"Bearer {token}",
        "X-GitHub-Api-Version": "2022-11-28",
    }


def get_lint_bot_comments(repo, token, pr_number):
    """Get the comments from the linting bot."""
    # repo is in the form of "org/repo"
    response = requests.get(
        f"https://api.github.com/repos/{repo}/issues/{pr_number}/comments",
        headers=get_headers(token),
    )
    response.raise_for_status()
    comments = response.json()

    failed_comment = "This PR is introducing linting issues. Here's a summary of the"
    success_comment = (
        "All linting checks passed. Your pull request is in excellent shape"
    )

    return [
        comment
        for comment in comments
        if comment["user"]["login"] == "github-actions[bot]"
        and (failed_comment in comment["body"] or success_comment in comment["body"])
    ]


def delete_existing_messages(comments, repo, token):
    """Delete the existing messages from the linting bot."""
    # repo is in the form of "org/repo"
    print("deleting comments")
    for comment in comments:
        response = requests.delete(
            f"https://api.github.com/repos/{repo}/issues/comments/{comment['id']}",
            headers=get_headers(token),
        )
        response.raise_for_status()


def create_comment(comment, repo, pr_number, token):
    """Create a new comment."""
    # repo is in the form of "org/repo"
    print("creating new comment")
    response = requests.post(
        f"https://api.github.com/repos/{repo}/issues/{pr_number}/comments",
        json={"body": comment},
        headers=get_headers(token),
    )
    response.raise_for_status()


if __name__ == "__main__":
    repo = os.environ["GITHUB_REPOSITORY"]
    token = os.environ["GITHUB_TOKEN"]
    pr_number = os.environ["PR_NUMBER"]

    if not repo or not token or not pr_number:
        raise ValueError(
            "One of the following environment variables is not set: "
            "GITHUB_REPOSITORY, GITHUB_TOKEN, PR_NUMBER"
        )

    delete_existing_messages(get_lint_bot_comments(repo, token, pr_number), repo, token)
    create_comment(message := get_message(), repo, pr_number, token)
    print(message)
