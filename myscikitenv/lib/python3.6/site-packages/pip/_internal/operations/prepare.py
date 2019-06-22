"""Prepares a distribution for installation
"""

import logging
import os

from pip._vendor import pkg_resources, requests

from pip._internal.build_env import BuildEnvironment
from pip._internal.download import (
    is_dir_url, is_file_url, is_vcs_url, unpack_url, url_to_path,
)
from pip._internal.exceptions import (
    DirectoryUrlHashUnsupported, HashUnpinned, InstallationError,
    PreviousBuildDirError, VcsHashUnsupported,
)
from pip._internal.utils.compat import expanduser
from pip._internal.utils.hashes import MissingHashes
from pip._internal.utils.logging import indent_log
from pip._internal.utils.misc import display_path, normalize_path
from pip._internal.utils.typing import MYPY_CHECK_RUNNING

if MYPY_CHECK_RUNNING:
    from typing import Any, Optional
    from pip._internal.req.req_install import InstallRequirement
    from pip._internal.index import PackageFinder
    from pip._internal.download import PipSession
    from pip._internal.req.req_tracker import RequirementTracker

logger = logging.getLogger(__name__)


def make_abstract_dist(req):
    # type: (InstallRequirement) -> DistAbstraction
    """Factory to make an abstract dist object.

    Preconditions: Either an editable req with a source_dir, or satisfied_by or
    a wheel link, or a non-editable req with a source_dir.

    :return: A concrete DistAbstraction.
    """
    if req.editable:
        return IsSDist(req)
    elif req.link and req.link.is_wheel:
        return IsWheel(req)
    else:
        return IsSDist(req)


class DistAbstraction(object):
    """Abstracts out the wheel vs non-wheel Resolver.resolve() logic.

    The requirements for anything installable are as follows:
     - we must be able to determine the requirement name
       (or we can't correctly handle the non-upgrade case).
     - we must be able to generate a list of run-time dependencies
       without installing any additional packages (or we would
       have to either burn time by doing temporary isolated installs
       or alternatively violate pips 'don't start installing unless
       all requirements are available' rule - neither of which are
       desirable).
     - for packages with setup requirements, we must also be able
       to determine their requirements without installing additional
       packages (for the same reason as run-time dependencies)
     - we must be able to create a Distribution object exposing the
       above metadata.
    """

    def __init__(self, req):
        # type: (InstallRequirement) -> None
        self.req = req  # type: InstallRequirement

    def dist(self):
        # type: () -> Any
        """Return a setuptools Dist object."""
        raise NotImplementedError

    def prep_for_dist(self, finder, build_isolation):
        # type: (PackageFinder, bool) -> Any
        """Ensure that we can get a Dist for this requirement."""
        raise NotImplementedError


class IsWheel(DistAbstraction):

    def dist(self):
        # type: () -> pkg_resources.Distribution
        return list(pkg_resources.find_distributions(
            self.req.source_dir))[0]

    def prep_for_dist(self, finder, build_isolation):
        # type: (PackageFinder, bool) -> Any
        # FIXME:https://github.com/pypa/pip/issues/1112
        pass


class IsSDist(DistAbstraction):

    def dist(self):
        return self.req.get_dist()

    def _raise_conflicts(self, conflicting_with, conflicting_reqs):
        conflict_messages = [
            '%s is incompatible with %s' % (installed, wanted)
            for installed, wanted in sorted(conflicting_reqs)
        ]
        raise InstallationError(
            "Some build dependencies for %s conflict with %s: %s." % (
                self.req, conflicting_with, ', '.join(conflict_messages))
        )

    def install_backend_dependencies(self, finder):
        # type: (PackageFinder) -> None
        """
        Install any extra build dependencies that the backend requests.

        :param finder: a PackageFinder object.
        """
        req = self.req
        with req.build_env:
            # We need to have the env active when calling the hook.
            req.spin_message = "Getting requirements to build wheel"
            reqs = req.pep517_backend.get_requires_for_build_wheel()
        conflicting, missing = req.build_env.check_requirements(reqs)
        if conflicting:
            self._raise_conflicts("the backend dependencies", conflicting)
        req.build_env.install_requirements(
            finder, missing, 'normal',
            "Installing backend dependencies"
        )

    def prep_for_dist(self, finder, build_isolation):
        # type: (PackageFinder, bool) -> None
        # Prepare for building. We need to:
        #   1. Load pyproject.toml (if it exists)
        #   2. Set up the build environment

        self.req.load_pyproject_toml()
        should_isolate = self.req.use_pep517 and build_isolation

        if should_isolate:
            # Isolate in a BuildEnvironment and install the build-time
            # requirements.
            self.req.build_env = BuildEnvironment()
            self.req.build_env.install_requirements(
                finder, self.req.pyproject_requires, 'overlay',
                "Installing build dependencies"
            )
            conflicting, missing = self.req.build_env.check_requirements(
                self.req.requirements_to_check
            )
            if conflicting:
                self._raise_conflicts("PEP 517/518 supported requirements",
                                      conflicting)
            if missing:
                logger.warning(
                    "Missing build requirements in pyproject.toml for %s.",
                    self.req,
                )
                logger.warning(
                    "The project does not specify a build backend, and "
                    "pip cannot fall back to setuptools without %s.",
                    " and ".join(map(repr, sorted(missing)))
                )

            # Install any extra build dependencies that the backend requests.
            # This must be done in a second pass, as the pyproject.toml
            # dependencies must be installed before we can call the backend.
            self.install_backend_dependencies(finder=finder)

        self.req.prepare_metadata()
        self.req.assert_source_matches_version()


class Installed(DistAbstraction):

    def dist(self):
        # type: () -> pkg_resources.Distribution
        return self.req.satisfied_by

    def prep_for_dist(self, finder, build_isolation):
        # type: (PackageFinder, bool) -> Any
        pass


class RequirementPreparer(object):
    """Prepares a Requirement
    """

    def __init__(
        self,
        build_dir,  # type: str
        download_dir,  # type: Optional[str]
        src_dir,  # type: str
        wheel_download_dir,  # type: Optional[str]
        progress_bar,  # type: str
        build_isolation,  # type: bool
        req_tracker  # type: RequirementTracker
    ):
        # type: (...) -> None
        super(RequirementPreparer, self).__init__()

        self.src_dir = src_dir
        self.build_dir = build_dir
        self.req_tracker = req_tracker

        # Where still packed archives should be written to. If None, they are
        # not saved, and are deleted immediately after unpacking.
        self.download_dir = download_dir

        # Where still-packed .whl files should be written to. If None, they are
        # written to the download_dir parameter. Separate to download_dir to
        # permit only keeping wheel archives for pip wheel.
        if wheel_download_dir:
            wheel_download_dir = normalize_path(wheel_download_dir)
        self.wheel_download_dir = wheel_download_dir

        # NOTE
        # download_dir and wheel_download_dir overlap semantically and may
        # be combined if we're willing to have non-wheel archives present in
        # the wheelhouse output by 'pip wheel'.

        self.progress_bar = progress_bar

        # Is build isolation allowed?
        self.build_isolation = build_isolation

    @property
    def _download_should_save(self):
        # type: () -> bool
        # TODO: Modify to reduce indentation needed
        if self.download_dir:
            self.download_dir = expanduser(self.download_dir)
            if os.path.exists(self.download_dir):
                return True
            else:
                logger.critical('Could not find download directory')
                raise InstallationError(
                    "Could not find or access download directory '%s'"
                    % display_path(self.download_dir))
        return False

    def prepare_linked_requirement(
        self,
        req,  # type: InstallRequirement
        session,  # type: PipSession
        finder,  # type: PackageFinder
        upgrade_allowed,  # type: bool
        require_hashes  # type: bool
    ):
        # type: (...) -> DistAbstraction
        """Prepare a requirement that would be obtained from req.link
        """
        # TODO: Breakup into smaller functions
        if req.link and req.link.scheme == 'file':
            path = url_to_path(req.link.url)
            logger.info('Processing %s', display_path(path))
        else:
            logger.info('Collecting %s', req)

        with indent_log():
            # @@ if filesystem packages are not marked
            # editable in a req, a non deterministic error
            # occurs when the script attempts to unpack the
            # build directory
            req.ensure_has_source_dir(self.build_dir)
            # If a checkout exists, it's unwise to keep going.  version
            # inconsistencies are logged later, but do not fail the
            # installation.
            # FIXME: this won't upgrade when there's an existing
            # package unpacked in `req.source_dir`
            # package unpacked in `req.source_dir`
            if os.path.exists(os.path.join(req.source_dir, 'setup.py')):
                raise PreviousBuildDirError(
                    "pip can't proceed with requirements '%s' due to a"
                    " pre-existing build directory (%s). This is "
                    "likely due to a previous installation that failed"
                    ". pip is being responsible and not assuming it "
                    "can delete this. Please delete it and try again."
                    % (req, req.source_dir)
                )
            req.populate_link(finder, upgrade_allowed, require_hashes)

            # We can't hit this spot and have populate_link return None.
            # req.satisfied_by is None here (because we're
            # guarded) and upgrade has no impact except when satisfied_by
            # is not None.
            # Then inside find_requirement existing_applicable -> False
            # If no new versions are found, DistributionNotFound is raised,
            # otherwise a result is guaranteed.
            assert req.link
            link = req.link

            # Now that we have the real link, we can tell what kind of
            # requirements we have and raise some more informative errors
            # than otherwise. (For example, we can raise VcsHashUnsupported
            # for a VCS URL rather than HashMissing.)
            if require_hashes:
                # We could check these first 2 conditions inside
                # unpack_url and save repetition of conditions, but then
                # we would report less-useful error messages for
                # unhashable requirements, complaining that there's no
                # hash provided.
                if is_vcs_url(link):
                    raise VcsHashUnsupported()
                elif is_file_url(link) and is_dir_url(link):
                    raise DirectoryUrlHashUnsupported()
                if not req.original_link and not req.is_pinned:
                    # Unpinned packages are asking for trouble when a new
                    # version is uploaded. This isn't a security check, but
                    # it saves users a surprising hash mismatch in the
                    # future.
                    #
                    # file:/// URLs aren't pinnable, so don't complain
                    # about them not being pinned.
                    raise HashUnpinned()

            hashes = req.hashes(trust_internet=not require_hashes)
            if require_hashes and not hashes:
                # Known-good hashes are missing for this requirement, so
                # shim it with a facade object that will provoke hash
                # computation and then raise a HashMissing exception
                # showing the user what the hash should be.
                hashes = MissingHashes()

            try:
                download_dir = self.download_dir
                # We always delete unpacked sdists after pip ran.
                autodelete_unpacked = True
                if req.link.is_wheel and self.wheel_download_dir:
                    # when doing 'pip wheel` we download wheels to a
                    # dedicated dir.
                    download_dir = self.wheel_download_dir
                if req.link.is_wheel:
                    if download_dir:
                        # When downloading, we only unpack wheels to get
                        # metadata.
                        autodelete_unpacked = True
                    else:
                        # When installing a wheel, we use the unpacked
                        # wheel.
                        autodelete_unpacked = False
                unpack_url(
                    req.link, req.source_dir,
                    download_dir, autodelete_unpacked,
                    session=session, hashes=hashes,
                    progress_bar=self.progress_bar
                )
            except requests.HTTPError as exc:
                logger.critical(
                    'Could not install requirement %s because of error %s',
                    req,
                    exc,
                )
                raise InstallationError(
                    'Could not install requirement %s because of HTTP '
                    'error %s for URL %s' %
                    (req, exc, req.link)
                )
            abstract_dist = make_abstract_dist(req)
            with self.req_tracker.track(req):
                abstract_dist.prep_for_dist(finder, self.build_isolation)
            if self._download_should_save:
                # Make a .zip of the source_dir we already created.
                if not req.link.is_artifact:
                    req.archive(self.download_dir)
        return abstract_dist

    def prepare_editable_requirement(
        self,
        req,  # type: InstallRequirement
        require_hashes,  # type: bool
        use_user_site,  # type: bool
        finder  # type: PackageFinder
    ):
        # type: (...) -> DistAbstraction
        """Prepare an editable requirement
        """
        assert req.editable, "cannot prepare a non-editable req as editable"

        logger.info('Obtaining %s', req)

        with indent_log():
            if require_hashes:
                raise InstallationError(
                    'The editable requirement %s cannot be installed when '
                    'requiring hashes, because there is no single file to '
                    'hash.' % req
                )
            req.ensure_has_source_dir(self.src_dir)
            req.update_editable(not self._download_should_save)

            abstract_dist = make_abstract_dist(req)
            with self.req_tracker.track(req):
                abstract_dist.prep_for_dist(finder, self.build_isolation)

            if self._download_should_save:
                req.archive(self.download_dir)
            req.check_if_exists(use_user_site)

        return abstract_dist

    def prepare_installed_requirement(self, req, require_hashes, skip_reason):
        # type: (InstallRequirement, bool, Optional[str]) -> DistAbstraction
        """Prepare an already-installed requirement
        """
        assert req.satisfied_by, "req should have been satisfied but isn't"
        assert skip_reason is not None, (
            "did not get skip reason skipped but req.satisfied_by "
            "is set to %r" % (req.satisfied_by,)
        )
        logger.info(
            'Requirement %s: %s (%s)',
            skip_reason, req, req.satisfied_by.version
        )
        with indent_log():
            if require_hashes:
                logger.debug(
                    'Since it is already installed, we are trusting this '
                    'package without checking its hash. To ensure a '
                    'completely repeatable environment, install into an '
                    'empty virtualenv.'
                )
            abstract_dist = Installed(req)

        return abstract_dist
