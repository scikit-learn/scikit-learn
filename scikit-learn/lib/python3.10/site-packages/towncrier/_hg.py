# Copyright (c) towncrier contributors, 2025
# See LICENSE for details.

from __future__ import annotations

import os

from collections.abc import Container
from subprocess import STDOUT, call, check_output


def get_default_compare_branch(branches: Container[str]) -> str | None:
    if "default" in branches:
        return "default"
    return None


def remove_files(fragment_filenames: list[str]) -> None:
    if not fragment_filenames:
        return

    # Filter out files that are unknown to mercurial
    hg_fragments = (
        check_output(["hg", "files"] + fragment_filenames, encoding="utf-8")
        .strip()
        .split("\n")
    )

    hg_fragments = [os.path.abspath(f) for f in hg_fragments if os.path.isfile(f)]
    fragment_filenames = [
        os.path.abspath(f) for f in fragment_filenames if os.path.isfile(f)
    ]
    call(["hg", "rm", "--force"] + hg_fragments, encoding="utf-8")
    unknown_fragments = set(fragment_filenames) - set(hg_fragments)
    for unknown_fragment in unknown_fragments:
        os.remove(unknown_fragment)


def stage_newsfile(directory: str, filename: str) -> None:
    call(["hg", "add", os.path.join(directory, filename)])


def get_remote_branches(base_directory: str) -> list[str]:
    branches = check_output(
        ["hg", "branches", "--template", "{branch}\n"],
        cwd=base_directory,
        encoding="utf-8",
    ).splitlines()

    return branches


def list_changed_files_compared_to_branch(
    base_directory: str, compare_with: str, include_staged: bool
) -> list[str]:
    output = check_output(
        ["hg", "diff", "--stat", "-r", compare_with],
        cwd=base_directory,
        encoding="utf-8",
        stderr=STDOUT,
    ).splitlines()

    return [line.split("|")[0].strip() for line in output if "|" in line]
