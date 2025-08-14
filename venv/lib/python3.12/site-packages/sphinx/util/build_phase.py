"""Build phase of Sphinx application."""

from __future__ import annotations

from enum import IntEnum


class BuildPhase(IntEnum):
    """Build phase of Sphinx application."""

    INITIALIZATION = 1
    READING = 2
    CONSISTENCY_CHECK = 3
    RESOLVING = 3
    WRITING = 4
