# Copyright (c) Stephen Finucane, 2019
# See LICENSE for details.

"""
Entry point of the command line interface.

Each sub-command has its separate CLI definition andd help messages.
"""

from __future__ import annotations

import click

from .build import _main as _build_cmd
from .check import _main as _check_cmd
from .click_default_group import DefaultGroup
from .create import _main as _create_cmd


@click.group(cls=DefaultGroup, default="build", default_if_no_args=True)
@click.version_option()
def cli() -> None:
    """
    Towncrier is a utility to produce useful, summarised news files for your project.
    Rather than reading the Git history as some newer tools to produce it, or having
    one single file which developers all write to, towncrier reads "news fragments"
    which contain information useful to end users.

    Towncrier delivers the news which is convenient to those that hear it, not those
    that write it.

    That is, a “news fragment” (a small file containing just enough information to
    be useful to end users) can be written that summarises what has changed from the
    “developer log” (which may contain complex information about the original issue,
    how it was fixed, who authored the fix, and who reviewed the fix). By compiling
    a collection of these fragments, towncrier can produce a digest of the changes
    which is valuable to those who may wish to use the software.
    """
    pass


cli.add_command(_build_cmd)
cli.add_command(_check_cmd)
cli.add_command(_create_cmd)
