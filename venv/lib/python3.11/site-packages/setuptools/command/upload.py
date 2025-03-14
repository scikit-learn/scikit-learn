from setuptools.dist import Distribution
from setuptools.errors import RemovedCommandError

from distutils import log
from distutils.command import upload as orig


class upload(orig.upload):
    """Formerly used to upload packages to PyPI."""

    distribution: Distribution  # override distutils.dist.Distribution with setuptools.dist.Distribution

    def run(self):
        msg = (
            "The upload command has been removed, use twine to upload "
            "instead (https://pypi.org/p/twine)"
        )

        self.announce("ERROR: " + msg, log.ERROR)
        raise RemovedCommandError(msg)
