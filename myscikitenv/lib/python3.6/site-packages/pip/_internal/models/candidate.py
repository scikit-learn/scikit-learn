from pip._vendor.packaging.version import parse as parse_version

from pip._internal.utils.models import KeyBasedCompareMixin
from pip._internal.utils.typing import MYPY_CHECK_RUNNING

if MYPY_CHECK_RUNNING:
    from pip._vendor.packaging.version import _BaseVersion
    from pip._internal.models.link import Link
    from typing import Any


class InstallationCandidate(KeyBasedCompareMixin):
    """Represents a potential "candidate" for installation.
    """

    def __init__(self, project, version, location):
        # type: (Any, str, Link) -> None
        self.project = project
        self.version = parse_version(version)  # type: _BaseVersion
        self.location = location

        super(InstallationCandidate, self).__init__(
            key=(self.project, self.version, self.location),
            defining_class=InstallationCandidate
        )

    def __repr__(self):
        # type: () -> str
        return "<InstallationCandidate({!r}, {!r}, {!r})>".format(
            self.project, self.version, self.location,
        )
