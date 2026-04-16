# SPDX-License-Identifier: MIT

"""
This module focuses on reading pyproject.toml fields with error collection. It is
mostly internal, except for License and Readme classes, which are re-exported in
the top-level package.

These functions return None if an error occurred, and use the ErrorCollector to
record errors.
"""

from __future__ import annotations

import dataclasses
import typing
from typing import Any

import packaging.requirements

if typing.TYPE_CHECKING:
    import pathlib
    from collections.abc import Generator, Iterable

    from packaging.requirements import Requirement

    from .errors import ErrorCollector
    from .project_table import ProjectTable


__all__ = [
    "License",
    "Readme",
]


def __dir__() -> list[str]:
    return __all__


@dataclasses.dataclass(frozen=True)
class License:
    """
    This represents a classic license, which contains text, and optionally a
    file path. Modern licenses are just SPDX identifiers, which are strings.
    """

    text: str
    file: pathlib.Path | None


@dataclasses.dataclass(frozen=True)
class Readme:
    """
    This represents a readme, which contains text and a content type, and
    optionally a file path.
    """

    text: str
    file: pathlib.Path | None
    content_type: str


def ensure_str(value: object) -> str | None:
    """Ensure that a value is a string."""
    if isinstance(value, str):
        return value
    return None


def ensure_list(val: object) -> list[str] | None:
    """Ensure that a value is a list of strings."""
    if val is None:
        return None
    if not isinstance(val, list):
        return None
    for item in val:
        if not isinstance(item, str):
            return None

    return val


def ensure_dict(val: object) -> dict[str, str] | None:
    """Ensure that a value is a dictionary of strings."""
    if not isinstance(val, dict):
        return None
    for item in val.values():
        if not isinstance(item, str):
            return None
    return val


def ensure_people(val: object) -> list[tuple[str, str | None]]:
    """Ensure that a value is a list of tables with optional "name" and "email" keys."""
    if not isinstance(val, list):
        return []
    for each in val:
        if not isinstance(each, dict):
            return []
        for value in each.values():
            if not isinstance(value, str):
                return []
    return [(entry.get("name", "Unknown"), entry.get("email")) for entry in val]


def get_license(
    project: ProjectTable, project_dir: pathlib.Path, error_collector: ErrorCollector
) -> License | str | None:
    """Get the license field from the project table. Handles PEP 639 style license too.

    None is returned if the license field is not present or if an error occurred.
    """
    val = project.get("license")
    if val is None:
        return None
    if isinstance(val, str):
        return val

    if isinstance(val, dict):
        _license = ensure_dict(val)
        if _license is None:
            return None
    else:
        return None

    file: pathlib.Path | None = None
    filename = _license.get("file")
    text = _license.get("text")

    if (filename and text) or (not filename and not text):
        return None

    if filename:
        file = project_dir.joinpath(filename)
        if not file.is_file():
            msg = f"License file not found ({filename!r})"
            error_collector.config_error(msg, key="project.license.file")
            return None
        text = file.read_text(encoding="utf-8")

    assert text is not None
    return License(text, file)


def get_license_files(
    project: dict[str, Any], project_dir: pathlib.Path, error_collector: ErrorCollector
) -> list[pathlib.Path] | None:
    """Get the license-files list of files from the project table.

    Returns None if an error occurred (including invalid globs, etc) or if
    not present.
    """
    license_files = project.get("license-files")
    if license_files is None:
        return None
    if ensure_list(license_files) is None:
        return None

    return list(_get_files_from_globs(project_dir, license_files, error_collector))


def get_readme(  # noqa: C901
    project: dict[str, Any], project_dir: pathlib.Path, error_collector: ErrorCollector
) -> Readme | None:
    """Get the text of the readme from the project table.

    Returns None if an error occurred or if the readme field is not present.
    """
    if "readme" not in project:
        return None

    filename: str | None = None
    file: pathlib.Path | None = None
    text: str | None = None
    content_type: str | None = None

    readme = project["readme"]
    if isinstance(readme, str):
        # readme is a file
        text = None
        filename = readme
        if filename.endswith(".md"):
            content_type = "text/markdown"
        elif filename.endswith(".rst"):
            content_type = "text/x-rst"
        else:
            msg = "Could not infer content type for readme file {filename!r}"
            error_collector.config_error(msg, key="project.readme", filename=filename)
            return None
    elif isinstance(readme, dict):
        # readme is a dict containing either 'file' or 'text', and content-type
        for field in readme:
            if field not in ("content-type", "file", "text"):
                return None

        content_type_raw = readme.get("content-type")
        if content_type_raw is not None:
            content_type = ensure_str(content_type_raw)
            if content_type is None:
                return None
        filename_raw = readme.get("file")
        if filename_raw is not None:
            filename = ensure_str(filename_raw)
            if filename is None:
                return None

        text_raw = readme.get("text")
        if text_raw is not None:
            text = ensure_str(text_raw)
            if text is None:
                return None

        if (filename and text) or (not filename and not text):
            return None
        if not content_type:
            return None
    else:
        return None

    if filename:
        file = project_dir.joinpath(filename)
        if not file.is_file():
            msg = "Readme file not found ({filename!r})"
            error_collector.config_error(
                msg, key="project.readme.file", filename=filename
            )
            return None
        text = file.read_text(encoding="utf-8")

    assert text is not None
    return Readme(text, file, content_type)


def get_dependencies(project: dict[str, Any]) -> list[Requirement]:
    """Get the dependencies from the project table."""
    requirement_strings: list[str] | None = None
    requirement_strings_raw = project.get("dependencies")
    if requirement_strings_raw is not None:
        requirement_strings = ensure_list(requirement_strings_raw)
    if requirement_strings is None:
        return []

    requirements: list[Requirement] = []
    try:
        requirements.extend(
            packaging.requirements.Requirement(req) for req in requirement_strings
        )
    except packaging.requirements.InvalidRequirement:
        return []
    return requirements


def get_optional_dependencies(
    project: dict[str, Any],
) -> dict[str, list[Requirement]]:
    """Get the optional dependencies from the project table."""
    val = project.get("optional-dependencies")
    if not val:
        return {}

    requirements_dict: dict[str, list[Requirement]] = {}
    if not isinstance(val, dict):
        return {}
    for extra, requirements in val.copy().items():
        assert isinstance(extra, str)
        if not isinstance(requirements, list):
            return {}
        requirements_dict[extra] = []
        for req in requirements:
            if not isinstance(req, str):
                return {}
            try:
                requirements_dict[extra].append(packaging.requirements.Requirement(req))
            except packaging.requirements.InvalidRequirement:
                return {}
    return dict(requirements_dict)


def get_entrypoints(project: dict[str, Any]) -> dict[str, dict[str, str]]:
    """Get the entrypoints from the project table."""
    val = project.get("entry-points")
    if val is None:
        return {}
    if not isinstance(val, dict):
        return {}
    for entrypoints in val.values():
        if not isinstance(entrypoints, dict):
            return {}
        for entrypoint in entrypoints.values():
            if not isinstance(entrypoint, str):
                return {}
    return val


def get_dynamic(project: dict[str, Any]) -> list[str]:
    """Get the dynamic fields from the project table.

    Returns an empty list if the field is not present or if an error occurred.
    """
    return ensure_list(project.get("dynamic", [])) or []


def _get_files_from_globs(
    project_dir: pathlib.Path, globs: Iterable[str], error_collector: ErrorCollector
) -> Generator[pathlib.Path, None, None]:
    """Given a list of globs, get files that match."""
    for glob in globs:
        if glob.startswith(("..", "/")):
            msg = "{glob!r} is an invalid {key} glob: the pattern must match files within the project directory"
            error_collector.config_error(msg, key="project.license-files", glob=glob)
            break
        files = [f for f in project_dir.glob(glob) if f.is_file()]
        if not files:
            msg = "Every pattern in {key} must match at least one file: {glob!r} did not match any"
            error_collector.config_error(msg, key="project.license-files", glob=glob)
            break
        for f in files:
            yield f.relative_to(project_dir)
