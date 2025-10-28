# SPDX-License-Identifier: MIT

"""
This is pyproject_metadata, a library for working with PEP 621 metadata.

Example usage:

.. code-block:: python

   from pyproject_metadata import StandardMetadata

   metadata = StandardMetadata.from_pyproject(
       parsed_pyproject, allow_extra_keys=False, all_errors=True, metadata_version="2.3"
   )

   pkg_info = metadata.as_rfc822()
   with open("METADATA", "wb") as f:
       f.write(pkg_info.as_bytes())

   ep = self.metadata.entrypoints.copy()
   ep["console_scripts"] = self.metadata.scripts
   ep["gui_scripts"] = self.metadata.gui_scripts
   for group, entries in ep.items():
       if entries:
           with open("entry_points.txt", "w", encoding="utf-8") as f:
               print(f"[{group}]", file=f)
               for name, target in entries.items():
                   print(f"{name} = {target}", file=f)
               print(file=f)

"""

from __future__ import annotations

import copy
import dataclasses
import email.message
import email.policy
import email.utils
import os
import os.path
import pathlib
import sys
import typing
import warnings

# Build backends may vendor this package, so all imports are relative.
from . import constants
from .errors import ConfigurationError, ConfigurationWarning, ErrorCollector
from .pyproject import License, PyProjectReader, Readme

if typing.TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any

    from packaging.requirements import Requirement

    if sys.version_info < (3, 11):
        from typing_extensions import Self
    else:
        from typing import Self

    from .project_table import Dynamic, PyProjectTable

import packaging.markers
import packaging.specifiers
import packaging.utils
import packaging.version

if sys.version_info < (3, 12, 4):
    import re

    RE_EOL_STR = re.compile(r"[\r\n]+")
    RE_EOL_BYTES = re.compile(rb"[\r\n]+")


__version__ = "0.9.1"

__all__ = [
    "ConfigurationError",
    "License",
    "RFC822Message",
    "RFC822Policy",
    "Readme",
    "StandardMetadata",
    "extras_build_system",
    "extras_project",
    "extras_top_level",
    "field_to_metadata",
]


def __dir__() -> list[str]:
    return __all__


def field_to_metadata(field: str) -> frozenset[str]:
    """
    Return the METADATA fields that correspond to a project field.
    """
    return frozenset(constants.PROJECT_TO_METADATA[field])


def extras_top_level(pyproject_table: Mapping[str, Any]) -> set[str]:
    """
    Return any extra keys in the top-level of the pyproject table.
    """
    return set(pyproject_table) - constants.KNOWN_TOPLEVEL_FIELDS


def extras_build_system(pyproject_table: Mapping[str, Any]) -> set[str]:
    """
    Return any extra keys in the build-system table.
    """
    return (
        set(pyproject_table.get("build-system", []))
        - constants.KNOWN_BUILD_SYSTEM_FIELDS
    )


def extras_project(pyproject_table: Mapping[str, Any]) -> set[str]:
    """
    Return any extra keys in the project table.
    """
    return set(pyproject_table.get("project", [])) - constants.KNOWN_PROJECT_FIELDS


@dataclasses.dataclass
class _SmartMessageSetter:
    """
    This provides a nice internal API for setting values in an Message to
    reduce boilerplate.

    If a value is None, do nothing.
    """

    message: email.message.Message

    def __setitem__(self, name: str, value: str | None) -> None:
        if not value:
            return
        self.message[name] = value

    def set_payload(self, payload: str) -> None:
        self.message.set_payload(payload)


@dataclasses.dataclass
class _JSonMessageSetter:
    """
    This provides an API to build a JSON message output in the same way as the
    classic Message. Line breaks are preserved this way.
    """

    data: dict[str, str | list[str]]

    def __setitem__(self, name: str, value: str | None) -> None:
        name = name.lower()
        key = name.replace("-", "_")

        if value is None:
            return

        if name == "keywords":
            values = (x.strip() for x in value.split(","))
            self.data[key] = [x for x in values if x]
        elif name in constants.KNOWN_MULTIUSE:
            entry = self.data.setdefault(key, [])
            assert isinstance(entry, list)
            entry.append(value)
        else:
            self.data[key] = value

    def set_payload(self, payload: str) -> None:
        self["description"] = payload


class RFC822Policy(email.policy.EmailPolicy):
    """
    This is :class:`email.policy.EmailPolicy`, but with a simple ``header_store_parse``
    implementation that handles multiline values, and some nice defaults.
    """

    utf8 = True
    mangle_from_ = False
    max_line_length = 0

    def header_store_parse(self, name: str, value: str) -> tuple[str, str]:
        if name.lower() not in constants.KNOWN_METADATA_FIELDS:
            msg = f"Unknown field {name!r}"
            raise ConfigurationError(msg, key=name)
        size = len(name) + 2
        value = value.replace("\n", "\n" + " " * size)
        return (name, value)

    if sys.version_info < (3, 12, 4):
        # Work around Python bug https://github.com/python/cpython/issues/117313
        def _fold(
            self, name: str, value: Any, refold_binary: bool = False
        ) -> str:  # pragma: no cover
            if hasattr(value, "name"):
                return value.fold(policy=self)  # type: ignore[no-any-return]
            maxlen = self.max_line_length if self.max_line_length else sys.maxsize

            # this is from the library version, and it improperly breaks on chars like 0x0c, treating
            # them as 'form feed' etc.
            # we need to ensure that only CR/LF is used as end of line
            # lines = value.splitlines()

            # this is a workaround which splits only on CR/LF characters
            if isinstance(value, bytes):
                lines = RE_EOL_BYTES.split(value)
            else:
                lines = RE_EOL_STR.split(value)

            refold = self.refold_source == "all" or (
                self.refold_source == "long"
                and (
                    (lines and len(lines[0]) + len(name) + 2 > maxlen)
                    or any(len(x) > maxlen for x in lines[1:])
                )
            )
            if refold or (refold_binary and email.policy._has_surrogates(value)):  # type: ignore[attr-defined]
                return self.header_factory(name, "".join(lines)).fold(policy=self)  # type: ignore[arg-type,no-any-return]
            return name + ": " + self.linesep.join(lines) + self.linesep  # type: ignore[arg-type]


class RFC822Message(email.message.EmailMessage):
    """
    This is :class:`email.message.EmailMessage` with two small changes: it defaults to
    our `RFC822Policy`, and it correctly writes unicode when being called
    with `bytes()`.
    """

    def __init__(self) -> None:
        super().__init__(policy=RFC822Policy())

    def as_bytes(
        self, unixfrom: bool = False, policy: email.policy.Policy | None = None
    ) -> bytes:
        """
        This handles unicode encoding.
        """
        return self.as_string(unixfrom, policy=policy).encode("utf-8")


@dataclasses.dataclass
class StandardMetadata:
    """
    This class represents the standard metadata fields for a project. It can be
    used to read metadata from a pyproject.toml table, validate it, and write it
    to an RFC822 message or JSON.
    """

    name: str
    version: packaging.version.Version | None = None
    description: str | None = None
    license: License | str | None = None
    license_files: list[pathlib.Path] | None = None
    readme: Readme | None = None
    requires_python: packaging.specifiers.SpecifierSet | None = None
    dependencies: list[Requirement] = dataclasses.field(default_factory=list)
    optional_dependencies: dict[str, list[Requirement]] = dataclasses.field(
        default_factory=dict
    )
    entrypoints: dict[str, dict[str, str]] = dataclasses.field(default_factory=dict)
    authors: list[tuple[str, str | None]] = dataclasses.field(default_factory=list)
    maintainers: list[tuple[str, str | None]] = dataclasses.field(default_factory=list)
    urls: dict[str, str] = dataclasses.field(default_factory=dict)
    classifiers: list[str] = dataclasses.field(default_factory=list)
    keywords: list[str] = dataclasses.field(default_factory=list)
    scripts: dict[str, str] = dataclasses.field(default_factory=dict)
    gui_scripts: dict[str, str] = dataclasses.field(default_factory=dict)
    dynamic: list[Dynamic] = dataclasses.field(default_factory=list)
    """
    This field is used to track dynamic fields. You can't set a field not in this list.
    """

    dynamic_metadata: list[str] = dataclasses.field(default_factory=list)
    """
    This is a list of METADATA fields that can change in between SDist and wheel. Requires metadata_version 2.2+.
    """
    metadata_version: str | None = None
    """
    This is the target metadata version. If None, it will be computed as a minimum based on the fields set.
    """
    all_errors: bool = False
    """
    If True, all errors will be collected and raised in an ExceptionGroup.
    """

    def __post_init__(self) -> None:
        self.validate()

    @property
    def auto_metadata_version(self) -> str:
        """
        This computes the metadata version based on the fields set in the object
        if ``metadata_version`` is None.
        """
        if self.metadata_version is not None:
            return self.metadata_version

        if isinstance(self.license, str) or self.license_files is not None:
            return "2.4"
        if self.dynamic_metadata:
            return "2.2"
        return "2.1"

    @property
    def canonical_name(self) -> str:
        """
        Return the canonical name of the project.
        """
        return packaging.utils.canonicalize_name(self.name)

    @classmethod
    def from_pyproject(  # noqa: C901
        cls,
        data: Mapping[str, Any],
        project_dir: str | os.PathLike[str] = os.path.curdir,
        metadata_version: str | None = None,
        dynamic_metadata: list[str] | None = None,
        *,
        allow_extra_keys: bool | None = None,
        all_errors: bool = False,
    ) -> Self:
        """
        Read metadata from a pyproject.toml table. This is the main method for
        creating an instance of this class. It also supports two additional
        fields: ``allow_extra_keys`` to control what happens when extra keys are
        present in the pyproject table, and ``all_errors``, to  raise all errors
        in an ExceptionGroup instead of raising the first one.
        """
        pyproject = PyProjectReader(collect_errors=all_errors)

        pyproject_table: PyProjectTable = data  # type: ignore[assignment]
        if "project" not in pyproject_table:
            msg = "Section {key} missing in pyproject.toml"
            pyproject.config_error(msg, key="project")
            pyproject.finalize("Failed to parse pyproject.toml")
            msg = "Unreachable code"  # pragma: no cover
            raise AssertionError(msg)  # pragma: no cover

        project = pyproject_table["project"]
        project_dir = pathlib.Path(project_dir)

        if not allow_extra_keys:
            extra_keys = extras_project(data)
            if extra_keys:
                extra_keys_str = ", ".join(sorted(f"{k!r}" for k in extra_keys))
                msg = "Extra keys present in {key}: {extra_keys}"
                pyproject.config_error(
                    msg,
                    key="project",
                    extra_keys=extra_keys_str,
                    warn=allow_extra_keys is None,
                )

        dynamic = pyproject.get_dynamic(project)

        for field in dynamic:
            if field in data["project"]:
                msg = 'Field {key} declared as dynamic in "project.dynamic" but is defined'
                pyproject.config_error(msg, key=f"project.{field}")

        raw_name = project.get("name")
        name = "UNKNOWN"
        if raw_name is None:
            msg = "Field {key} missing"
            pyproject.config_error(msg, key="project.name")
        else:
            tmp_name = pyproject.ensure_str(raw_name, "project.name")
            if tmp_name is not None:
                name = tmp_name

        version: packaging.version.Version | None = packaging.version.Version("0.0.0")
        raw_version = project.get("version")
        if raw_version is not None:
            version_string = pyproject.ensure_str(raw_version, "project.version")
            if version_string is not None:
                try:
                    version = (
                        packaging.version.Version(version_string)
                        if version_string
                        else None
                    )
                except packaging.version.InvalidVersion:
                    msg = "Invalid {key} value, expecting a valid PEP 440 version"
                    pyproject.config_error(
                        msg, key="project.version", got=version_string
                    )
        elif "version" not in dynamic:
            msg = (
                "Field {key} missing and 'version' not specified in \"project.dynamic\""
            )
            pyproject.config_error(msg, key="project.version")

        # Description fills Summary, which cannot be multiline
        # However, throwing an error isn't backward compatible,
        # so leave it up to the users for now.
        project_description_raw = project.get("description")
        description = (
            pyproject.ensure_str(project_description_raw, "project.description")
            if project_description_raw is not None
            else None
        )

        requires_python_raw = project.get("requires-python")
        requires_python = None
        if requires_python_raw is not None:
            requires_python_string = pyproject.ensure_str(
                requires_python_raw, "project.requires-python"
            )
            if requires_python_string is not None:
                try:
                    requires_python = packaging.specifiers.SpecifierSet(
                        requires_python_string
                    )
                except packaging.specifiers.InvalidSpecifier:
                    msg = "Invalid {key} value, expecting a valid specifier set"
                    pyproject.config_error(
                        msg, key="project.requires-python", got=requires_python_string
                    )

        self = None
        with pyproject.collect():
            self = cls(
                name=name,
                version=version,
                description=description,
                license=pyproject.get_license(project, project_dir),
                license_files=pyproject.get_license_files(project, project_dir),
                readme=pyproject.get_readme(project, project_dir),
                requires_python=requires_python,
                dependencies=pyproject.get_dependencies(project),
                optional_dependencies=pyproject.get_optional_dependencies(project),
                entrypoints=pyproject.get_entrypoints(project),
                authors=pyproject.ensure_people(
                    project.get("authors", []), "project.authors"
                ),
                maintainers=pyproject.ensure_people(
                    project.get("maintainers", []), "project.maintainers"
                ),
                urls=pyproject.ensure_dict(project.get("urls", {}), "project.urls")
                or {},
                classifiers=pyproject.ensure_list(
                    project.get("classifiers", []), "project.classifiers"
                )
                or [],
                keywords=pyproject.ensure_list(
                    project.get("keywords", []), "project.keywords"
                )
                or [],
                scripts=pyproject.ensure_dict(
                    project.get("scripts", {}), "project.scripts"
                )
                or {},
                gui_scripts=pyproject.ensure_dict(
                    project.get("gui-scripts", {}), "project.gui-scripts"
                )
                or {},
                dynamic=dynamic,
                dynamic_metadata=dynamic_metadata or [],
                metadata_version=metadata_version,
                all_errors=all_errors,
            )

        pyproject.finalize("Failed to parse pyproject.toml")
        assert self is not None
        return self

    def as_rfc822(self) -> RFC822Message:
        """
        Return an RFC822 message with the metadata.
        """
        message = RFC822Message()
        smart_message = _SmartMessageSetter(message)
        self._write_metadata(smart_message)
        return message

    def as_json(self) -> dict[str, str | list[str]]:
        """
        Return a JSON message with the metadata.
        """
        message: dict[str, str | list[str]] = {}
        smart_message = _JSonMessageSetter(message)
        self._write_metadata(smart_message)
        return message

    def validate(self, *, warn: bool = True) -> None:  # noqa: C901
        """
        Validate metadata for consistency and correctness. Will also produce
        warnings if ``warn`` is given. Respects ``all_errors``. This is called
        when loading a pyproject.toml, and when making metadata. Checks:

        - ``metadata_version`` is a known version or None
        - ``name`` is a valid project name
        - ``license_files`` can't be used with classic ``license``
        - License classifiers can't be used with SPDX license
        - ``description`` is a single line (warning)
        - ``license`` is not an SPDX license expression if metadata_version >= 2.4 (warning)
        - License classifiers deprecated for metadata_version >= 2.4 (warning)
        - ``license`` is an SPDX license expression if metadata_version >= 2.4
        - ``license_files`` is supported only for metadata_version >= 2.4
        - ``project_url`` can't contain keys over 32 characters
        """
        errors = ErrorCollector(collect_errors=self.all_errors)

        if self.auto_metadata_version not in constants.KNOWN_METADATA_VERSIONS:
            msg = "The metadata_version must be one of {versions} or None (default)"
            errors.config_error(msg, versions=constants.KNOWN_METADATA_VERSIONS)

        try:
            packaging.utils.canonicalize_name(self.name, validate=True)
        except packaging.utils.InvalidName:
            msg = (
                "Invalid project name {name!r}. A valid name consists only of ASCII letters and "
                "numbers, period, underscore and hyphen. It must start and end with a letter or number"
            )
            errors.config_error(msg, key="project.name", name=self.name)

        if self.license_files is not None and isinstance(self.license, License):
            msg = '{key} must not be used when "project.license" is not a SPDX license expression'
            errors.config_error(msg, key="project.license-files")

        if isinstance(self.license, str) and any(
            c.startswith("License ::") for c in self.classifiers
        ):
            msg = "Setting {key} to an SPDX license expression is not compatible with 'License ::' classifiers"
            errors.config_error(msg, key="project.license")

        if warn:
            if self.description and "\n" in self.description:
                warnings.warn(
                    'The one-line summary "project.description" should not contain more than one line. Readers might merge or truncate newlines.',
                    ConfigurationWarning,
                    stacklevel=2,
                )
            if self.auto_metadata_version not in constants.PRE_SPDX_METADATA_VERSIONS:
                if isinstance(self.license, License):
                    warnings.warn(
                        'Set "project.license" to an SPDX license expression for metadata >= 2.4',
                        ConfigurationWarning,
                        stacklevel=2,
                    )
                elif any(c.startswith("License ::") for c in self.classifiers):
                    warnings.warn(
                        "'License ::' classifiers are deprecated for metadata >= 2.4, use a SPDX license expression for \"project.license\" instead",
                        ConfigurationWarning,
                        stacklevel=2,
                    )

        if (
            isinstance(self.license, str)
            and self.auto_metadata_version in constants.PRE_SPDX_METADATA_VERSIONS
        ):
            msg = "Setting {key} to an SPDX license expression is supported only when emitting metadata version >= 2.4"
            errors.config_error(msg, key="project.license")

        if (
            self.license_files is not None
            and self.auto_metadata_version in constants.PRE_SPDX_METADATA_VERSIONS
        ):
            msg = "{key} is supported only when emitting metadata version >= 2.4"
            errors.config_error(msg, key="project.license-files")

        for name in self.urls:
            if len(name) > 32:
                msg = "{key} names cannot be more than 32 characters long"
                errors.config_error(msg, key="project.urls", got=name)

        errors.finalize("Metadata validation failed")

    def _write_metadata(  # noqa: C901
        self, smart_message: _SmartMessageSetter | _JSonMessageSetter
    ) -> None:
        """
        Write the metadata to the message. Handles JSON or Message.
        """
        errors = ErrorCollector(collect_errors=self.all_errors)
        with errors.collect():
            self.validate(warn=False)

        smart_message["Metadata-Version"] = self.auto_metadata_version
        smart_message["Name"] = self.name
        if not self.version:
            msg = "Field {key} missing"
            errors.config_error(msg, key="project.version")
        smart_message["Version"] = str(self.version)
        # skip 'Platform'
        # skip 'Supported-Platform'
        if self.description:
            smart_message["Summary"] = self.description
        smart_message["Keywords"] = ",".join(self.keywords) or None
        # skip 'Home-page'
        # skip 'Download-URL'
        smart_message["Author"] = _name_list(self.authors)
        smart_message["Author-Email"] = _email_list(self.authors)
        smart_message["Maintainer"] = _name_list(self.maintainers)
        smart_message["Maintainer-Email"] = _email_list(self.maintainers)

        if isinstance(self.license, License):
            smart_message["License"] = self.license.text
        elif isinstance(self.license, str):
            smart_message["License-Expression"] = self.license

        if self.license_files is not None:
            for license_file in sorted(set(self.license_files)):
                smart_message["License-File"] = os.fspath(license_file.as_posix())
        elif (
            self.auto_metadata_version not in constants.PRE_SPDX_METADATA_VERSIONS
            and isinstance(self.license, License)
            and self.license.file
        ):
            smart_message["License-File"] = os.fspath(self.license.file.as_posix())

        for classifier in self.classifiers:
            smart_message["Classifier"] = classifier
        # skip 'Provides-Dist'
        # skip 'Obsoletes-Dist'
        # skip 'Requires-External'
        for name, url in self.urls.items():
            smart_message["Project-URL"] = f"{name}, {url}"
        if self.requires_python:
            smart_message["Requires-Python"] = str(self.requires_python)
        for dep in self.dependencies:
            smart_message["Requires-Dist"] = str(dep)
        for extra, requirements in self.optional_dependencies.items():
            norm_extra = extra.replace(".", "-").replace("_", "-").lower()
            smart_message["Provides-Extra"] = norm_extra
            for requirement in requirements:
                smart_message["Requires-Dist"] = str(
                    _build_extra_req(norm_extra, requirement)
                )
        if self.readme:
            if self.readme.content_type:
                smart_message["Description-Content-Type"] = self.readme.content_type
            smart_message.set_payload(self.readme.text)
        # Core Metadata 2.2
        if self.auto_metadata_version != "2.1":
            for field in self.dynamic_metadata:
                if field.lower() in {"name", "version", "dynamic"}:
                    msg = f"Metadata field {field!r} cannot be declared dynamic"
                    errors.config_error(msg)
                if field.lower() not in constants.KNOWN_METADATA_FIELDS:
                    msg = f"Unknown metadata field {field!r} cannot be declared dynamic"
                    errors.config_error(msg)
                smart_message["Dynamic"] = field

        errors.finalize("Failed to write metadata")


def _name_list(people: list[tuple[str, str | None]]) -> str | None:
    """
    Build a comma-separated list of names.
    """
    return ", ".join(name for name, email_ in people if not email_) or None


def _email_list(people: list[tuple[str, str | None]]) -> str | None:
    """
    Build a comma-separated list of emails.
    """
    return (
        ", ".join(
            email.utils.formataddr((name, _email)) for name, _email in people if _email
        )
        or None
    )


def _build_extra_req(
    extra: str,
    requirement: Requirement,
) -> Requirement:
    """
    Build a new requirement with an extra marker.
    """
    requirement = copy.copy(requirement)
    if requirement.marker:
        if "or" in requirement.marker._markers:
            requirement.marker = packaging.markers.Marker(
                f"({requirement.marker}) and extra == {extra!r}"
            )
        else:
            requirement.marker = packaging.markers.Marker(
                f"{requirement.marker} and extra == {extra!r}"
            )
    else:
        requirement.marker = packaging.markers.Marker(f"extra == {extra!r}")
    return requirement
