from setuptools.errors import RemovedCommandError

from ..dist import Distribution

import distutils.command.register as orig
from distutils import log


class register(orig.register):
    """Formerly used to register packages on PyPI."""

    distribution: Distribution  # override distutils.dist.Distribution with setuptools.dist.Distribution

    def run(self):
        msg = (
            "The register command has been removed, use twine to upload "
            "instead (https://pypi.org/p/twine)"
        )

        self.announce("ERROR: " + msg, log.ERROR)

        raise RemovedCommandError(msg)
