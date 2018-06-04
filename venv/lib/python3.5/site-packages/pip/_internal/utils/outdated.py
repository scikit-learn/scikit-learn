from __future__ import absolute_import

import datetime
import json
import logging
import os.path
import sys

from pip._vendor import lockfile, pkg_resources
from pip._vendor.packaging import version as packaging_version

from pip._internal.compat import WINDOWS
from pip._internal.index import PackageFinder
from pip._internal.utils.filesystem import check_path_owner
from pip._internal.utils.misc import ensure_dir, get_installed_version

SELFCHECK_DATE_FMT = "%Y-%m-%dT%H:%M:%SZ"


logger = logging.getLogger(__name__)


class SelfCheckState(object):
    def __init__(self, cache_dir):
        self.statefile_path = os.path.join(cache_dir, "selfcheck.json")

        # Load the existing state
        try:
            with open(self.statefile_path) as statefile:
                self.state = json.load(statefile)[sys.prefix]
        except (IOError, ValueError, KeyError):
            self.state = {}

    def save(self, pypi_version, current_time):
        # Check to make sure that we own the directory
        if not check_path_owner(os.path.dirname(self.statefile_path)):
            return

        # Now that we've ensured the directory is owned by this user, we'll go
        # ahead and make sure that all our directories are created.
        ensure_dir(os.path.dirname(self.statefile_path))

        # Attempt to write out our version check file
        with lockfile.LockFile(self.statefile_path):
            if os.path.exists(self.statefile_path):
                with open(self.statefile_path) as statefile:
                    state = json.load(statefile)
            else:
                state = {}

            state[sys.prefix] = {
                "last_check": current_time.strftime(SELFCHECK_DATE_FMT),
                "pypi_version": pypi_version,
            }

            with open(self.statefile_path, "w") as statefile:
                json.dump(state, statefile, sort_keys=True,
                          separators=(",", ":"))


def was_installed_by_pip(pkg):
    """Checks whether pkg was installed by pip

    This is used not to display the upgrade message when pip is in fact
    installed by system package manager, such as dnf on Fedora.
    """
    try:
        dist = pkg_resources.get_distribution(pkg)
        return (dist.has_metadata('INSTALLER') and
                'pip' in dist.get_metadata_lines('INSTALLER'))
    except pkg_resources.DistributionNotFound:
        return False


def pip_version_check(session, options):
    """Check for an update for pip.

    Limit the frequency of checks to once per week. State is stored either in
    the active virtualenv or in the user's USER_CACHE_DIR keyed off the prefix
    of the pip script path.
    """
    installed_version = get_installed_version("pip")
    if not installed_version:
        return

    pip_version = packaging_version.parse(installed_version)
    pypi_version = None

    try:
        state = SelfCheckState(cache_dir=options.cache_dir)

        current_time = datetime.datetime.utcnow()
        # Determine if we need to refresh the state
        if "last_check" in state.state and "pypi_version" in state.state:
            last_check = datetime.datetime.strptime(
                state.state["last_check"],
                SELFCHECK_DATE_FMT
            )
            if (current_time - last_check).total_seconds() < 7 * 24 * 60 * 60:
                pypi_version = state.state["pypi_version"]

        # Refresh the version if we need to or just see if we need to warn
        if pypi_version is None:
            # Lets use PackageFinder to see what the latest pip version is
            finder = PackageFinder(
                find_links=options.find_links,
                index_urls=[options.index_url] + options.extra_index_urls,
                allow_all_prereleases=False,  # Explicitly set to False
                trusted_hosts=options.trusted_hosts,
                process_dependency_links=options.process_dependency_links,
                session=session,
            )
            all_candidates = finder.find_all_candidates("pip")
            if not all_candidates:
                return
            pypi_version = str(
                max(all_candidates, key=lambda c: c.version).version
            )

            # save that we've performed a check
            state.save(pypi_version, current_time)

        remote_version = packaging_version.parse(pypi_version)

        # Determine if our pypi_version is older
        if (pip_version < remote_version and
                pip_version.base_version != remote_version.base_version and
                was_installed_by_pip('pip')):
            # Advise "python -m pip" on Windows to avoid issues
            # with overwriting pip.exe.
            if WINDOWS:
                pip_cmd = "python -m pip"
            else:
                pip_cmd = "pip"
            logger.warning(
                "You are using pip version %s, however version %s is "
                "available.\nYou should consider upgrading via the "
                "'%s install --upgrade pip' command.",
                pip_version, pypi_version, pip_cmd
            )
    except Exception:
        logger.debug(
            "There was an error checking the latest version of pip",
            exc_info=True,
        )
