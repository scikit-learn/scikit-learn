# SPDX-License-Identifier: MIT

"""
This module focues on reading pyproject.toml fields with error collection. It is
mostly internal, except for License and Readme classes, which are re-exported in
the top-level package.
"""

from __future__ import annotations

import dataclasses
import pathlib
import re
import typing

import packaging.requirements

from .errors import ErrorCollector

if typing.TYPE_CHECKING:
    from collections.abc import Generator, Iterable, Sequence

    from packaging.requirements import Requirement

    from .project_table import ContactTable, Dynamic, ProjectTable


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


T = typing.TypeVar("T")


@dataclasses.dataclass
class PyProjectReader(ErrorCollector):
    """Class for reading pyproject.toml fields with error collection.

    Unrelated errors are collected and raised at once if the `collect_errors`
    parameter is set to `True`. Some methods will return None if an error was
    raised. Most of them expect a non-None value as input to enforce the caller
    to handle missing vs. error correctly. The exact design is based on usage,
    as this is an internal class.
    """

    def ensure_str(self, value: str, key: str) -> str | None:
        """Ensure that a value is a string."""
        if isinstance(value, str):
            return value

        msg = "Field {key} has an invalid type, expecting a string"
        self.config_error(msg, key=key, got_type=type(value))
        return None

    def ensure_list(self, val: list[T], key: str) -> list[T] | None:
        """Ensure that a value is a list of strings."""
        if not isinstance(val, list):
            msg = "Field {key} has an invalid type, expecting a list of strings"
            self.config_error(msg, key=key, got_type=type(val))
            return None
        for item in val:
            if not isinstance(item, str):
                msg = "Field {key} contains item with invalid type, expecting a string"
                self.config_error(msg, key=key, got_type=type(item))
                return None

        return val

    def ensure_dict(self, val: dict[str, str], key: str) -> dict[str, str] | None:
        """Ensure that a value is a dictionary of strings."""
        if not isinstance(val, dict):
            msg = "Field {key} has an invalid type, expecting a table of strings"
            self.config_error(msg, key=key, got_type=type(val))
            return None
        for subkey, item in val.items():
            if not isinstance(item, str):
                msg = "Field {key} has an invalid type, expecting a string"
                self.config_error(msg, key=f"{key}.{subkey}", got_type=type(item))
                return None
        return val

    def ensure_people(
        self, val: Sequence[ContactTable], key: str
    ) -> list[tuple[str, str | None]]:
        """Ensure that a value is a list of tables with optional "name" and "email" keys."""
        if not isinstance(val, list):
            msg = (
                "Field {key} has an invalid type, expecting a list of "
                'tables containing the "name" and/or "email" keys'
            )
            self.config_error(msg, key=key, got_type=type(val))
            return []
        for each in val:
            if not isinstance(each, dict):
                msg = (
                    "Field {key} has an invalid type, expecting a list of "
                    'tables containing the "name" and/or "email" keys'
                    " (got list with {type_name})"
                )
                self.config_error(msg, key=key, type_name=type(each).__name__)
                return []
            for value in each.values():
                if not isinstance(value, str):
                    msg = (
                        "Field {key} has an invalid type, expecting a list of "
                        'tables containing the "name" and/or "email" keys'
                        " (got list with dict with {type_name})"
                    )
                    self.config_error(msg, key=key, type_name=type(value).__name__)
                    return []
            extra_keys = set(each) - {"name", "email"}
            if extra_keys:
                msg = (
                    "Field {key} has an invalid type, expecting a list of "
                    'tables containing the "name" and/or "email" keys'
                    " (got list with dict with extra keys {extra_keys})"
                )
                self.config_error(
                    msg,
                    key=key,
                    extra_keys=", ".join(sorted(f'"{k}"' for k in extra_keys)),
                )
                return []
        return [(entry.get("name", "Unknown"), entry.get("email")) for entry in val]

    def get_license(
        self, project: ProjectTable, project_dir: pathlib.Path
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
            _license = self.ensure_dict(val, "project.license")  # type: ignore[arg-type]
            if _license is None:
                return None
        else:
            msg = "Field {key} has an invalid type, expecting a string or table of strings"
            self.config_error(msg, key="project.license", got_type=type(val))
            return None

        for field in _license:
            if field not in ("file", "text"):
                msg = "Unexpected field {key}"
                self.config_error(msg, key=f"project.license.{field}")
                return None

        file: pathlib.Path | None = None
        filename = _license.get("file")
        text = _license.get("text")

        if (filename and text) or (not filename and not text):
            msg = (
                'Invalid {key} contents, expecting a string or one key "file" or "text"'
            )
            self.config_error(msg, key="project.license", got=_license)
            return None

        if filename:
            file = project_dir.joinpath(filename)
            if not file.is_file():
                msg = f"License file not found ({filename!r})"
                self.config_error(msg, key="project.license.file")
                return None
            text = file.read_text(encoding="utf-8")

        assert text is not None
        return License(text, file)

    def get_license_files(
        self, project: ProjectTable, project_dir: pathlib.Path
    ) -> list[pathlib.Path] | None:
        """Get the license-files list of files from the project table.

        Returns None if an error occurred (including invalid globs, etc) or if
        not present.
        """
        license_files = project.get("license-files")
        if license_files is None:
            return None
        if self.ensure_list(license_files, "project.license-files") is None:
            return None

        return list(self._get_files_from_globs(project_dir, license_files))

    def get_readme(  # noqa: C901
        self, project: ProjectTable, project_dir: pathlib.Path
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
                self.config_error(msg, key="project.readme", filename=filename)
                return None
        elif isinstance(readme, dict):
            # readme is a dict containing either 'file' or 'text', and content-type
            for field in readme:
                if field not in ("content-type", "file", "text"):
                    msg = "Unexpected field {key}"
                    self.config_error(msg, key=f"project.readme.{field}")
                    return None

            content_type_raw = readme.get("content-type")
            if content_type_raw is not None:
                content_type = self.ensure_str(
                    content_type_raw, "project.readme.content-type"
                )
                if content_type is None:
                    return None
            filename_raw = readme.get("file")
            if filename_raw is not None:
                filename = self.ensure_str(filename_raw, "project.readme.file")
                if filename is None:
                    return None

            text_raw = readme.get("text")
            if text_raw is not None:
                text = self.ensure_str(text_raw, "project.readme.text")
                if text is None:
                    return None

            if (filename and text) or (not filename and not text):
                msg = 'Invalid {key} contents, expecting either "file" or "text"'
                self.config_error(msg, key="project.readme", got=readme)
                return None
            if not content_type:
                msg = "Field {key} missing"
                self.config_error(msg, key="project.readme.content-type")
                return None
        else:
            msg = "Field {key} has an invalid type, expecting either a string or table of strings"
            self.config_error(msg, key="project.readme", got_type=type(readme))
            return None

        if filename:
            file = project_dir.joinpath(filename)
            if not file.is_file():
                msg = "Readme file not found ({filename!r})"
                self.config_error(msg, key="project.readme.file", filename=filename)
                return None
            text = file.read_text(encoding="utf-8")

        assert text is not None
        return Readme(text, file, content_type)

    def get_dependencies(self, project: ProjectTable) -> list[Requirement]:
        """Get the dependencies from the project table."""

        requirement_strings: list[str] | None = None
        requirement_strings_raw = project.get("dependencies")
        if requirement_strings_raw is not None:
            requirement_strings = self.ensure_list(
                requirement_strings_raw, "project.dependencies"
            )
        if requirement_strings is None:
            return []

        requirements: list[Requirement] = []
        for req in requirement_strings:
            try:
                requirements.append(packaging.requirements.Requirement(req))
            except packaging.requirements.InvalidRequirement as e:
                msg = "Field {key} contains an invalid PEP 508 requirement string {req!r} ({error!r})"
                self.config_error(msg, key="project.dependencies", req=req, error=e)
                return []
        return requirements

    def get_optional_dependencies(
        self,
        project: ProjectTable,
    ) -> dict[str, list[Requirement]]:
        """Get the optional dependencies from the project table."""

        val = project.get("optional-dependencies")
        if not val:
            return {}

        requirements_dict: dict[str, list[Requirement]] = {}
        if not isinstance(val, dict):
            msg = "Field {key} has an invalid type, expecting a table of PEP 508 requirement strings"
            self.config_error(
                msg, key="project.optional-dependencies", got_type=type(val)
            )
            return {}
        for extra, requirements in val.copy().items():
            assert isinstance(extra, str)
            if not isinstance(requirements, list):
                msg = "Field {key} has an invalid type, expecting a table of PEP 508 requirement strings"
                self.config_error(
                    msg,
                    key=f"project.optional-dependencies.{extra}",
                    got_type=type(requirements),
                )
                return {}
            requirements_dict[extra] = []
            for req in requirements:
                if not isinstance(req, str):
                    msg = "Field {key} has an invalid type, expecting a PEP 508 requirement string"
                    self.config_error(
                        msg,
                        key=f"project.optional-dependencies.{extra}",
                        got_type=type(req),
                    )
                    return {}
                try:
                    requirements_dict[extra].append(
                        packaging.requirements.Requirement(req)
                    )
                except packaging.requirements.InvalidRequirement as e:
                    msg = (
                        "Field {key} contains "
                        "an invalid PEP 508 requirement string {req!r} ({error!r})"
                    )
                    self.config_error(
                        msg,
                        key=f"project.optional-dependencies.{extra}",
                        req=req,
                        error=e,
                    )
                    return {}
        return dict(requirements_dict)

    def get_entrypoints(self, project: ProjectTable) -> dict[str, dict[str, str]]:
        """Get the entrypoints from the project table."""

        val = project.get("entry-points", None)
        if val is None:
            return {}
        if not isinstance(val, dict):
            msg = "Field {key} has an invalid type, expecting a table of entrypoint sections"
            self.config_error(msg, key="project.entry-points", got_type=type(val))
            return {}
        for section, entrypoints in val.items():
            assert isinstance(section, str)
            if not re.match(r"^\w+(\.\w+)*$", section):
                msg = (
                    "Field {key} has an invalid value, expecting a name "
                    "containing only alphanumeric, underscore, or dot characters"
                )
                self.config_error(msg, key="project.entry-points", got=section)
                return {}
            if not isinstance(entrypoints, dict):
                msg = (
                    "Field {key} has an invalid type, expecting a table of entrypoints"
                )
                self.config_error(
                    msg,
                    key=f"project.entry-points.{section}",
                    got_type=type(entrypoints),
                )
                return {}
            for name, entrypoint in entrypoints.items():
                assert isinstance(name, str)
                if not isinstance(entrypoint, str):
                    msg = "Field {key} has an invalid type, expecting a string"
                    self.config_error(
                        msg,
                        key=f"project.entry-points.{section}.{name}",
                        got_type=type(entrypoint),
                    )
                    return {}
        return val

    def get_dynamic(self, project: ProjectTable) -> list[Dynamic]:
        """Get the dynamic fields from the project table.

        Returns an empty list if the field is not present or if an error occurred.
        """
        dynamic = project.get("dynamic", [])

        self.ensure_list(dynamic, "project.dynamic")

        if "name" in dynamic:
            msg = "Unsupported field 'name' in {key}"
            self.config_error(msg, key="project.dynamic")
            return []

        return dynamic

    def _get_files_from_globs(
        self, project_dir: pathlib.Path, globs: Iterable[str]
    ) -> Generator[pathlib.Path, None, None]:
        """Given a list of globs, get files that match."""

        for glob in globs:
            if glob.startswith(("..", "/")):
                msg = "{glob!r} is an invalid {key} glob: the pattern must match files within the project directory"
                self.config_error(msg, key="project.license-files", glob=glob)
                break
            files = [f for f in project_dir.glob(glob) if f.is_file()]
            if not files:
                msg = "Every pattern in {key} must match at least one file: {glob!r} did not match any"
                self.config_error(msg, key="project.license-files", glob=glob)
                break
            for f in files:
                yield f.relative_to(project_dir)
