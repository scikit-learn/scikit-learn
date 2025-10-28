# Copyright (c) Amber Brown, 2018
# See LICENSE for details.


from __future__ import annotations

import os
import sys

from subprocess import CalledProcessError

import click

from ._builder import find_fragments
from ._settings import config_option_help, load_config_from_options
from ._vcs import (
    get_default_compare_branch,
    get_remote_branches,
    list_changed_files_compared_to_branch,
)


@click.command(name="check")
@click.option(
    "--compare-with",
    default=None,
    metavar="BRANCH",
    help=(
        "Checks files changed running git diff --name-only BRANCH... "
        "BRANCH is the branch to be compared with. "
        "Default to origin/main"
    ),
)
@click.option(
    "--dir",
    "directory",
    default=None,
    metavar="PATH",
    help="Check fragment in directory. Default to current directory.",
)
@click.option(
    "--config",
    "config",
    default=None,
    metavar="FILE_PATH",
    help=config_option_help,
)
@click.option(
    "--staged",
    "staged",
    is_flag=True,
    default=False,
    help="Include staged files as part of the branch checked in the --compare-with",
)
def _main(
    compare_with: str | None, directory: str | None, config: str | None, staged: bool
) -> None:
    """
    Check for new fragments on a branch.
    """
    __main(compare_with, directory, config, staged)


def __main(
    comparewith: str | None,
    directory: str | None,
    config_path: str | None,
    staged: bool,
) -> None:
    base_directory, config = load_config_from_options(directory, config_path)

    if comparewith is None:
        comparewith = get_default_compare_branch(
            base_directory, get_remote_branches(base_directory=base_directory)
        )

    if comparewith is None:
        click.echo("Could not detect default branch. Aborting.")
        sys.exit(1)

    try:
        files_changed = list_changed_files_compared_to_branch(
            base_directory, comparewith, staged
        )
    except CalledProcessError as e:
        click.echo("git produced output while failing:")
        click.echo(e.output)
        raise

    if not files_changed:
        click.echo(
            f"On {comparewith} branch, or no diffs, so no newsfragment required."
        )
        sys.exit(0)

    files = {os.path.abspath(path) for path in files_changed}

    click.echo("Looking at these files:")
    click.echo("----")
    for n, change in enumerate(files, start=1):
        click.echo(f"{n}. {change}")
    click.echo("----")

    # This will fail if any fragment files have an invalid name:
    all_fragment_files = find_fragments(base_directory, config, strict=True)[1]

    news_file = os.path.normpath(os.path.join(base_directory, config.filename))
    if news_file in files:
        click.echo("Checks SKIPPED: news file changes detected.")
        sys.exit(0)

    fragments = set()  # will only include fragments of types that are checked
    unchecked_fragments = set()  # will include fragments of types that are not checked
    for fragment_filename, category in all_fragment_files:
        if config.types[category]["check"]:
            fragments.add(fragment_filename)
        else:
            unchecked_fragments.add(fragment_filename)
    fragments_in_branch = fragments & files

    if not fragments_in_branch:
        unchecked_fragments_in_branch = unchecked_fragments & files
        if unchecked_fragments:
            click.echo("Found newsfragments of unchecked types in the branch:")
            for n, fragment in enumerate(unchecked_fragments_in_branch, start=1):
                click.echo(f"{n}. {fragment}")
        else:
            click.echo("No new newsfragments found on this branch.")
        sys.exit(1)
    else:
        click.echo("Found:")
        for n, fragment in enumerate(fragments_in_branch, start=1):
            click.echo(f"{n}. {fragment}")
        sys.exit(0)


if __name__ == "__main__":  # pragma: no cover
    _main()
