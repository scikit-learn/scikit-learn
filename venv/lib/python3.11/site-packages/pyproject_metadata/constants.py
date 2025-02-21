# SPDX-License-Identifier: MIT

"""
Constants for the pyproject_metadata package, collected here to make them easy
to update. These should be considered mostly private.
"""

from __future__ import annotations

__all__ = [
    "KNOWN_BUILD_SYSTEM_FIELDS",
    "KNOWN_METADATA_FIELDS",
    "KNOWN_METADATA_VERSIONS",
    "KNOWN_METADATA_VERSIONS",
    "KNOWN_MULTIUSE",
    "KNOWN_PROJECT_FIELDS",
    "KNOWN_TOPLEVEL_FIELDS",
    "PRE_SPDX_METADATA_VERSIONS",
    "PROJECT_TO_METADATA",
]


def __dir__() -> list[str]:
    return __all__


KNOWN_METADATA_VERSIONS = {"2.1", "2.2", "2.3", "2.4"}
PRE_SPDX_METADATA_VERSIONS = {"2.1", "2.2", "2.3"}

PROJECT_TO_METADATA = {
    "authors": frozenset(["Author", "Author-Email"]),
    "classifiers": frozenset(["Classifier"]),
    "dependencies": frozenset(["Requires-Dist"]),
    "description": frozenset(["Summary"]),
    "dynamic": frozenset(),
    "entry-points": frozenset(),
    "gui-scripts": frozenset(),
    "keywords": frozenset(["Keywords"]),
    "license": frozenset(["License", "License-Expression"]),
    "license-files": frozenset(["License-File"]),
    "maintainers": frozenset(["Maintainer", "Maintainer-Email"]),
    "name": frozenset(["Name"]),
    "optional-dependencies": frozenset(["Provides-Extra", "Requires-Dist"]),
    "readme": frozenset(["Description", "Description-Content-Type"]),
    "requires-python": frozenset(["Requires-Python"]),
    "scripts": frozenset(),
    "urls": frozenset(["Project-URL"]),
    "version": frozenset(["Version"]),
}

KNOWN_TOPLEVEL_FIELDS = {"build-system", "project", "tool", "dependency-groups"}
KNOWN_BUILD_SYSTEM_FIELDS = {"backend-path", "build-backend", "requires"}
KNOWN_PROJECT_FIELDS = set(PROJECT_TO_METADATA)

KNOWN_METADATA_FIELDS = {
    "author",
    "author-email",
    "classifier",
    "description",
    "description-content-type",
    "download-url",  # Not specified via pyproject standards, deprecated by PEP 753
    "dynamic",  # Can't be in dynamic
    "home-page",  # Not specified via pyproject standards, deprecated by PEP 753
    "keywords",
    "license",
    "license-expression",
    "license-file",
    "maintainer",
    "maintainer-email",
    "metadata-version",
    "name",  # Can't be in dynamic
    "obsoletes",  # Deprecated
    "obsoletes-dist",  # Rarely used
    "platform",  # Not specified via pyproject standards
    "project-url",
    "provides",  # Deprecated
    "provides-dist",  # Rarely used
    "provides-extra",
    "requires",  # Deprecated
    "requires-dist",
    "requires-external",  # Not specified via pyproject standards
    "requires-python",
    "summary",
    "supported-platform",  # Not specified via pyproject standards
    "version",  # Can't be in dynamic
}

KNOWN_MULTIUSE = {
    "dynamic",
    "platform",
    "provides-extra",
    "supported-platform",
    "license-file",
    "classifier",
    "requires-dist",
    "requires-external",
    "project-url",
    "provides-dist",
    "obsoletes-dist",
    "requires",  # Deprecated
    "obsoletes",  # Deprecated
    "provides",  # Deprecated
}
