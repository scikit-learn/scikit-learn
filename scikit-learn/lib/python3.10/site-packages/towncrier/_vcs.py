# Copyright (c) towncrier contributors, 2025
# See LICENSE for details.

from __future__ import annotations

import os

from collections.abc import Container
from typing import Protocol


class VCSMod(Protocol):
    def get_default_compare_branch(self, branches: Container[str]) -> str | None: ...
    def remove_files(self, fragment_filenames: list[str]) -> None: ...
    def stage_newsfile(self, directory: str, filename: str) -> None: ...
    def get_remote_branches(self, base_directory: str) -> list[str]: ...

    def list_changed_files_compared_to_branch(
        self, base_directory: str, compare_with: str, include_staged: bool
    ) -> list[str]: ...


def _get_mod(base_directory: str) -> VCSMod:
    base_directory = os.path.abspath(base_directory)
    if os.path.exists(os.path.join(base_directory, ".git")):
        from . import _git

        return _git
    elif os.path.exists(os.path.join(base_directory, ".hg")):
        from . import _hg

        hg: VCSMod = _hg

        return hg
    else:
        # No VCS was found in the current directory
        # We will try our luck in the parent directory.
        parent = os.path.dirname(base_directory)
        if parent == base_directory:
            # We reached the fs root, abandoning
            from . import _novcs

            return _novcs

        return _get_mod(parent)


def get_default_compare_branch(
    base_directory: str, branches: Container[str]
) -> str | None:
    return _get_mod(base_directory).get_default_compare_branch(branches)


def remove_files(base_directory: str, fragment_filenames: list[str]) -> None:
    return _get_mod(base_directory).remove_files(fragment_filenames)


def stage_newsfile(directory: str, filename: str) -> None:
    return _get_mod(directory).stage_newsfile(directory, filename)


def get_remote_branches(base_directory: str) -> list[str]:
    return _get_mod(base_directory).get_remote_branches(base_directory)


def list_changed_files_compared_to_branch(
    base_directory: str, compare_with: str, include_staged: bool
) -> list[str]:
    return _get_mod(base_directory).list_changed_files_compared_to_branch(
        base_directory,
        compare_with,
        include_staged,
    )
