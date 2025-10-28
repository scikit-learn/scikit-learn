# Copyright (c) Amber Brown, 2015
# See LICENSE for details.

from __future__ import annotations

import os

from collections.abc import Container
from subprocess import STDOUT, call, check_output
from warnings import warn


def get_default_compare_branch(branches: Container[str]) -> str | None:
    if "origin/main" in branches:
        return "origin/main"
    if "origin/master" in branches:
        warn(
            'Using "origin/master" as default compare branch is deprecated '
            "and will be removed in a future version.",
            DeprecationWarning,
            stacklevel=2,
        )
        return "origin/master"
    return None


def remove_files(fragment_filenames: list[str]) -> None:
    if not fragment_filenames:
        return

    # Filter out files that are unknown to git
    git_fragments = check_output(
        ["git", "ls-files"] + fragment_filenames, encoding="utf-8"
    ).split("\n")

    git_fragments = [os.path.abspath(f) for f in git_fragments if os.path.isfile(f)]
    call(["git", "rm", "--quiet", "--force"] + git_fragments)
    unknown_fragments = set(fragment_filenames) - set(git_fragments)
    for unknown_fragment in unknown_fragments:
        os.remove(unknown_fragment)


def stage_newsfile(directory: str, filename: str) -> None:
    call(["git", "add", os.path.join(directory, filename)])


def get_remote_branches(base_directory: str) -> list[str]:
    output = check_output(
        ["git", "branch", "-r"], cwd=base_directory, encoding="utf-8", stderr=STDOUT
    )

    return [branch.strip() for branch in output.strip().splitlines()]


def list_changed_files_compared_to_branch(
    base_directory: str, compare_with: str, include_staged: bool
) -> list[str]:
    output = check_output(
        ["git", "diff", "--name-only", compare_with + "..."],
        cwd=base_directory,
        encoding="utf-8",
        stderr=STDOUT,
    )
    filenames = output.strip().splitlines()
    if include_staged:
        output = check_output(
            ["git", "diff", "--name-only", "--cached"],
            cwd=base_directory,
            encoding="utf-8",
            stderr=STDOUT,
        )
        filenames.extend(output.strip().splitlines())

    return filenames
