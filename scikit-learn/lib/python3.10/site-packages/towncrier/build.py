# Copyright (c) Amber Brown, 2015
# See LICENSE for details.

"""
Build a combined news file from news fragments.
"""

from __future__ import annotations

import os
import re
import sys

from datetime import date
from pathlib import Path

import click

from click import Context, Option, UsageError

from towncrier import _vcs

from ._builder import find_fragments, render_fragments, split_fragments
from ._project import get_project_name, get_version
from ._settings import ConfigError, config_option_help, load_config_from_options
from ._writer import append_to_newsfile


if sys.version_info < (3, 10):
    import importlib_resources as resources
else:
    from importlib import resources


def _get_date() -> str:
    return date.today().isoformat()


def _validate_answer(ctx: Context, param: Option, value: bool) -> bool:
    value_check = (
        ctx.params.get("answer_yes")
        if param.name == "answer_keep"
        else ctx.params.get("answer_keep")
    )
    if value_check and value:
        click.echo("You can not choose both --yes and --keep at the same time")
        ctx.abort()
    return value


@click.command(name="build")
@click.option(
    "--draft",
    "draft",
    default=False,
    flag_value=True,
    help=(
        "Render the news fragments to standard output. "
        "Don't write to files, don't check versions."
    ),
)
@click.option(
    "--config",
    "config_file",
    default=None,
    metavar="FILE_PATH",
    help=config_option_help,
)
@click.option(
    "--dir",
    "directory",
    default=None,
    metavar="PATH",
    help="Build fragment in directory. Default to current directory.",
)
@click.option(
    "--name",
    "project_name",
    default=None,
    help="Pass a custom project name.",
)
@click.option(
    "--version",
    "project_version",
    default=None,
    help="Render the news fragments using given version.",
)
@click.option(
    "--date",
    "project_date",
    default=None,
    help="Render the news fragments using the given date.",
)
@click.option(
    "--yes",
    "answer_yes",
    default=None,
    flag_value=True,
    help="Do not ask for confirmation to remove news fragments.",
    callback=_validate_answer,
)
@click.option(
    "--keep",
    "answer_keep",
    default=None,
    flag_value=True,
    help="Do not ask for confirmations. But keep news fragments.",
    callback=_validate_answer,
)
def _main(
    draft: bool,
    directory: str | None,
    config_file: str | None,
    project_name: str | None,
    project_version: str | None,
    project_date: str | None,
    answer_yes: bool,
    answer_keep: bool,
) -> None:
    """
    Build a combined news file from news fragment.
    """
    try:
        return __main(
            draft,
            directory,
            config_file,
            project_name,
            project_version,
            project_date,
            answer_yes,
            answer_keep,
        )
    except ConfigError as e:
        print(e, file=sys.stderr)
        sys.exit(1)


def __main(
    draft: bool,
    directory: str | None,
    config_file: str | None,
    project_name: str | None,
    project_version: str | None,
    project_date: str | None,
    answer_yes: bool,
    answer_keep: bool,
) -> None:
    """
    The main entry point.
    """
    base_directory, config = load_config_from_options(directory, config_file)
    to_err = draft

    if project_version is None:
        project_version = config.version
    if project_version is None:
        if not config.package:
            raise UsageError(
                "'--version' is required since the config file does "
                "not contain 'version' or 'package'."
            )
        project_version = get_version(
            os.path.join(base_directory, config.package_dir), config.package
        ).strip()

    click.echo("Loading template...", err=to_err)
    if isinstance(config.template, tuple):
        template = (
            resources.files(config.template[0])
            .joinpath(config.template[1])
            .read_text(encoding="utf-8")
        )
        template_extension = os.path.splitext(config.template[1])[1]
    else:
        template = Path(config.template).read_text(encoding="utf-8")
        template_extension = os.path.splitext(config.template)[1]
    is_markdown = template_extension.lower() == ".md"

    click.echo("Finding news fragments...", err=to_err)

    fragment_contents, fragment_files = find_fragments(
        base_directory,
        config,
        # Fail if any fragment filenames are invalid only if ignore list is set
        # (this maintains backward compatibility):
        strict=(config.ignore is not None),
    )
    fragment_filenames = [filename for (filename, _category) in fragment_files]

    click.echo("Rendering news fragments...", err=to_err)
    fragments = split_fragments(
        fragment_contents, config.types, all_bullets=config.all_bullets
    )

    if project_name is None:
        project_name = config.name
        if not project_name:
            package = config.package
            if package:
                project_name = get_project_name(
                    os.path.abspath(os.path.join(base_directory, config.package_dir)),
                    package,
                )
            else:
                # Can't determine a project_name, but maybe it is not needed.
                project_name = ""

    if project_date is None:
        project_date = _get_date().strip()

    # Render the title in the template if the title format is set to "". It can
    # alternatively be set to False or a string, in either case it shouldn't be rendered
    # in the template.
    render_title = config.title_format == ""

    # Add format-specific context to the template
    md_header_level = 1
    if is_markdown:
        if config.title_format:
            m = re.search(r"^#+(?=\s)", config.title_format, re.MULTILINE)
            lvl = len(m[0]) if m else 0
        else:  # TODO: derive from template or make configurable?
            lvl = 1 if render_title else 0
        md_header_level = lvl

    rendered = render_fragments(
        # The 0th underline is used for the top line
        template,
        config.issue_format,
        fragments,
        config.types,
        config.underlines[1:],
        config.wrap,
        {"name": project_name, "version": project_version, "date": project_date},
        top_underline=config.underlines[0],
        all_bullets=config.all_bullets,
        render_title=render_title,
        md_header_level=md_header_level,
    )

    if config.title_format:
        top_line = config.title_format.format(
            name=project_name, version=project_version, project_date=project_date
        )
        if is_markdown:
            parts = [top_line]
        else:
            parts = [top_line, config.underlines[0] * len(top_line)]
        parts.append(rendered)
        content = "\n".join(parts)
    else:
        top_line = ""
        content = rendered

    if draft:
        click.echo(
            "Draft only -- nothing has been written.\n"
            "What is seen below is what would be written.\n",
            err=to_err,
        )
        click.echo(content)
        return

    click.echo("Writing to newsfile...", err=to_err)
    news_file = config.filename

    if config.single_file is False:
        # The release notes for each version are stored in a separate file.
        # The name of that file is generated based on the current version and project.
        news_file = news_file.format(
            name=project_name, version=project_version, project_date=project_date
        )

    append_to_newsfile(
        base_directory,
        news_file,
        config.start_string,
        top_line,
        content,
        single_file=config.single_file,
    )

    click.echo("Staging newsfile...", err=to_err)
    _vcs.stage_newsfile(base_directory, news_file)

    if should_remove_fragment_files(
        fragment_filenames,
        answer_yes,
        answer_keep,
    ):
        click.echo("Removing news fragments...", err=to_err)
        _vcs.remove_files(base_directory, fragment_filenames)

    click.echo("Done!", err=to_err)


def should_remove_fragment_files(
    fragment_filenames: list[str],
    answer_yes: bool,
    answer_keep: bool,
) -> bool:
    if not fragment_filenames:
        click.echo("No news fragments to remove. Skipping!")
        return False
    try:
        if answer_keep:
            click.echo("Keeping the following files:")
            # Not proceeding with the removal of the files.
            return False

        if answer_yes:
            click.echo("Removing the following files:")
        else:
            click.echo("I want to remove the following files:")
    finally:
        # Will always be printed, even for answer_keep to help with possible troubleshooting
        for filename in fragment_filenames:
            click.echo(filename)

    if answer_yes or click.confirm("Is it okay if I remove those files?", default=True):
        return True
    return False


if __name__ == "__main__":  # pragma: no cover
    _main()
