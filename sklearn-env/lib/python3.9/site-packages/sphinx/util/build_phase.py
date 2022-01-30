"""
    sphinx.util.build_phase
    ~~~~~~~~~~~~~~~~~~~~~~~

    Build phase of Sphinx application.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from enum import IntEnum


class BuildPhase(IntEnum):
    """Build phase of Sphinx application."""
    INITIALIZATION = 1
    READING = 2
    CONSISTENCY_CHECK = 3
    RESOLVING = 3
    WRITING = 4
