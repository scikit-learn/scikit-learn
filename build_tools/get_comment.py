# This script is used to generate a comment for a PR when linting issues are
# detected. It is used by the `Comment on failed linting` GitHub Action.

import os
import re

from github import Auth, Github, GithubException


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


def get_message(log_file, repo_str, pr_number, sha, run_id, details, versions):
    with open(log_file, "r") as f:
        log = f.read()

    sub_text = (
        "\n\n<sub> _Generated for commit:"
        f" [{sha[:7]}](https://github.com/{repo_str}/pull/{pr_number}/commits/{sha}). "
        "Link to the linter CI: [here]"
        f"(https://github.com/{repo_str}/actions/runs/{run_id})_ </sub>"
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
        # no issues detected, the linting succeeded
        return None

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
        + "https://scikit-learn.org/dev/developers/development_setup.html#set-up-pre-commit)"
        + ".\n\n"
        + "You can see the details of the linting issues under the `lint` job [here]"
        + f"(https://github.com/{repo_str}/actions/runs/{run_id})\n\n"
        + message
        + sub_text
    )

    return message


def find_lint_bot_comments(issue):
    """Get the comment from the linting bot."""

    failed_comment = "❌ Linting issues"

    for comment in issue.get_comments():
        if comment.user.login == "github-actions[bot]":
            if failed_comment in comment.body:
                return comment

    return None


def create_or_update_comment(comment, message, issue):
    """Create a new comment or update the existing linting comment."""

    if comment is not None:
        print("Updating existing comment")
        comment.edit(message)
    else:
        print("Creating new comment")
        issue.create_comment(message)


def update_linter_fails_label(linting_failed, issue):
    """Add or remove the label indicating that the linting has failed."""

    label = "CI:Linter failure"

    if linting_failed:
        issue.add_to_labels(label)

    else:
        try:
            issue.remove_from_labels(label)
        except GithubException as exception:
            # The exception is ignored if raised because the issue did not have the
            # label already
            if not exception.message == "Label does not exist":
                raise


if __name__ == "__main__":
    repo_str = os.environ["GITHUB_REPOSITORY"]
    token = os.environ["GITHUB_TOKEN"]
    pr_number = os.environ["PR_NUMBER"]
    sha = os.environ["BRANCH_SHA"]
    log_file = os.environ["LOG_FILE"]
    run_id = os.environ["RUN_ID"]
    versions_file = os.environ["VERSIONS_FILE"]

    versions = get_versions(versions_file)

    for var, val in [
        ("GITHUB_REPOSITORY", repo_str),
        ("GITHUB_TOKEN", token),
        ("PR_NUMBER", pr_number),
        ("LOG_FILE", log_file),
        ("RUN_ID", run_id),
    ]:
        if not val:
            raise ValueError(f"The following environment variable is not set: {var}")

    if not re.match(r"\d+$", pr_number):
        raise ValueError(f"PR_NUMBER should be a number, got {pr_number!r} instead")
    pr_number = int(pr_number)

    gh = Github(auth=Auth.Token(token))
    repo = gh.get_repo(repo_str)
    issue = repo.get_issue(number=pr_number)

    message = get_message(
        log_file,
        repo_str=repo_str,
        pr_number=pr_number,
        sha=sha,
        run_id=run_id,
        details=True,
        versions=versions,
    )

    update_linter_fails_label(
        linting_failed=message is not None,
        issue=issue,
    )

    comment = find_lint_bot_comments(issue)

    if message is None:  # linting succeeded
        if comment is not None:
            print("Deleting existing comment.")
            comment.delete()
    else:
        try:
            create_or_update_comment(comment, message, issue)
            print(message)
        except GithubException:
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
            create_or_update_comment(comment, message, issue)
            print(message)
