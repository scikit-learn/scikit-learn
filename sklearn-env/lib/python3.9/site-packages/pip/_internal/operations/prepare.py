"""Prepares a distribution for installation
"""

# The following comment should be removed at some point in the future.
# mypy: strict-optional=False

import logging
import mimetypes
import os
import shutil

from pip._vendor.packaging.utils import canonicalize_name
from pip._vendor.six import PY2

from pip._internal.distributions import make_distribution_for_install_requirement
from pip._internal.distributions.installed import InstalledDistribution
from pip._internal.exceptions import (
    DirectoryUrlHashUnsupported,
    HashMismatch,
    HashUnpinned,
    InstallationError,
    NetworkConnectionError,
    PreviousBuildDirError,
    VcsHashUnsupported,
)
from pip._internal.models.wheel import Wheel
from pip._internal.network.download import BatchDownloader, Downloader
from pip._internal.network.lazy_wheel import (
    HTTPRangeRequestUnsupported,
    dist_from_wheel_url,
)
from pip._internal.utils.filesystem import copy2_fixed
from pip._internal.utils.hashes import MissingHashes
from pip._internal.utils.logging import indent_log
from pip._internal.utils.misc import display_path, hide_url, path_to_display, rmtree
from pip._internal.utils.temp_dir import TempDirectory
from pip._internal.utils.typing import MYPY_CHECK_RUNNING
from pip._internal.utils.unpacking import unpack_file
from pip._internal.vcs import vcs

if MYPY_CHECK_RUNNING:
    from typing import Callable, Dict, Iterable, List, Optional, Tuple

    from mypy_extensions import TypedDict
    from pip._vendor.pkg_resources import Distribution

    from pip._internal.index.package_finder import PackageFinder
    from pip._internal.models.link import Link
    from pip._internal.network.session import PipSession
    from pip._internal.req.req_install import InstallRequirement
    from pip._internal.req.req_tracker import RequirementTracker
    from pip._internal.utils.hashes import Hashes

    if PY2:
        CopytreeKwargs = TypedDict(
            'CopytreeKwargs',
            {
                'ignore': Callable[[str, List[str]], List[str]],
                'symlinks': bool,
            },
            total=False,
        )
    else:
        CopytreeKwargs = TypedDict(
            'CopytreeKwargs',
            {
                'copy_function': Callable[[str, str], None],
                'ignore': Callable[[str, List[str]], List[str]],
                'ignore_dangling_symlinks': bool,
                'symlinks': bool,
            },
            total=False,
        )

logger = logging.getLogger(__name__)


def _get_prepared_distribution(
    req,  # type: InstallRequirement
    req_tracker,  # type: RequirementTracker
    finder,  # type: PackageFinder
    build_isolation,  # type: bool
):
    # type: (...) -> Distribution
    """Prepare a distribution for installation."""
    abstract_dist = make_distribution_for_install_requirement(req)
    with req_tracker.track(req):
        abstract_dist.prepare_distribution_metadata(finder, build_isolation)
    return abstract_dist.get_pkg_resources_distribution()


def unpack_vcs_link(link, location):
    # type: (Link, str) -> None
    vcs_backend = vcs.get_backend_for_scheme(link.scheme)
    assert vcs_backend is not None
    vcs_backend.unpack(location, url=hide_url(link.url))


class File(object):

    def __init__(self, path, content_type):
        # type: (str, Optional[str]) -> None
        self.path = path
        if content_type is None:
            self.content_type = mimetypes.guess_type(path)[0]
        else:
            self.content_type = content_type


def get_http_url(
    link,  # type: Link
    download,  # type: Downloader
    download_dir=None,  # type: Optional[str]
    hashes=None,  # type: Optional[Hashes]
):
    # type: (...) -> File
    temp_dir = TempDirectory(kind="unpack", globally_managed=True)
    # If a download dir is specified, is the file already downloaded there?
    already_downloaded_path = None
    if download_dir:
        already_downloaded_path = _check_download_dir(
            link, download_dir, hashes
        )

    if already_downloaded_path:
        from_path = already_downloaded_path
        content_type = None
    else:
        # let's download to a tmp dir
        from_path, content_type = download(link, temp_dir.path)
        if hashes:
            hashes.check_against_path(from_path)

    return File(from_path, content_type)


def _copy2_ignoring_special_files(src, dest):
    # type: (str, str) -> None
    """Copying special files is not supported, but as a convenience to users
    we skip errors copying them. This supports tools that may create e.g.
    socket files in the project source directory.
    """
    try:
        copy2_fixed(src, dest)
    except shutil.SpecialFileError as e:
        # SpecialFileError may be raised due to either the source or
        # destination. If the destination was the cause then we would actually
        # care, but since the destination directory is deleted prior to
        # copy we ignore all of them assuming it is caused by the source.
        logger.warning(
            "Ignoring special file error '%s' encountered copying %s to %s.",
            str(e),
            path_to_display(src),
            path_to_display(dest),
        )


def _copy_source_tree(source, target):
    # type: (str, str) -> None
    target_abspath = os.path.abspath(target)
    target_basename = os.path.basename(target_abspath)
    target_dirname = os.path.dirname(target_abspath)

    def ignore(d, names):
        # type: (str, List[str]) -> List[str]
        skipped = []  # type: List[str]
        if d == source:
            # Pulling in those directories can potentially be very slow,
            # exclude the following directories if they appear in the top
            # level dir (and only it).
            # See discussion at https://github.com/pypa/pip/pull/6770
            skipped += ['.tox', '.nox']
        if os.path.abspath(d) == target_dirname:
            # Prevent an infinite recursion if the target is in source.
            # This can happen when TMPDIR is set to ${PWD}/...
            # and we copy PWD to TMPDIR.
            skipped += [target_basename]
        return skipped

    kwargs = dict(ignore=ignore, symlinks=True)  # type: CopytreeKwargs

    if not PY2:
        # Python 2 does not support copy_function, so we only ignore
        # errors on special file copy in Python 3.
        kwargs['copy_function'] = _copy2_ignoring_special_files

    shutil.copytree(source, target, **kwargs)


def get_file_url(
    link,  # type: Link
    download_dir=None,  # type: Optional[str]
    hashes=None  # type: Optional[Hashes]
):
    # type: (...) -> File
    """Get file and optionally check its hash.
    """
    # If a download dir is specified, is the file already there and valid?
    already_downloaded_path = None
    if download_dir:
        already_downloaded_path = _check_download_dir(
            link, download_dir, hashes
        )

    if already_downloaded_path:
        from_path = already_downloaded_path
    else:
        from_path = link.file_path

    # If --require-hashes is off, `hashes` is either empty, the
    # link's embedded hash, or MissingHashes; it is required to
    # match. If --require-hashes is on, we are satisfied by any
    # hash in `hashes` matching: a URL-based or an option-based
    # one; no internet-sourced hash will be in `hashes`.
    if hashes:
        hashes.check_against_path(from_path)
    return File(from_path, None)


def unpack_url(
    link,  # type: Link
    location,  # type: str
    download,  # type: Downloader
    download_dir=None,  # type: Optional[str]
    hashes=None,  # type: Optional[Hashes]
):
    # type: (...) -> Optional[File]
    """Unpack link into location, downloading if required.

    :param hashes: A Hashes object, one of whose embedded hashes must match,
        or HashMismatch will be raised. If the Hashes is empty, no matches are
        required, and unhashable types of requirements (like VCS ones, which
        would ordinarily raise HashUnsupported) are allowed.
    """
    # non-editable vcs urls
    if link.is_vcs:
        unpack_vcs_link(link, location)
        return None

    # If it's a url to a local directory
    if link.is_existing_dir():
        if os.path.isdir(location):
            rmtree(location)
        _copy_source_tree(link.file_path, location)
        return None

    # file urls
    if link.is_file:
        file = get_file_url(link, download_dir, hashes=hashes)

    # http urls
    else:
        file = get_http_url(
            link,
            download,
            download_dir,
            hashes=hashes,
        )

    # unpack the archive to the build dir location. even when only downloading
    # archives, they have to be unpacked to parse dependencies, except wheels
    if not link.is_wheel:
        unpack_file(file.path, location, file.content_type)

    return file


def _check_download_dir(link, download_dir, hashes):
    # type: (Link, str, Optional[Hashes]) -> Optional[str]
    """ Check download_dir for previously downloaded file with correct hash
        If a correct file is found return its path else None
    """
    download_path = os.path.join(download_dir, link.filename)

    if not os.path.exists(download_path):
        return None

    # If already downloaded, does its hash match?
    logger.info('File was already downloaded %s', download_path)
    if hashes:
        try:
            hashes.check_against_path(download_path)
        except HashMismatch:
            logger.warning(
                'Previously-downloaded file %s has bad hash. '
                'Re-downloading.',
                download_path
            )
            os.unlink(download_path)
            return None
    return download_path


class RequirementPreparer(object):
    """Prepares a Requirement
    """

    def __init__(
        self,
        build_dir,  # type: str
        download_dir,  # type: Optional[str]
        src_dir,  # type: str
        build_isolation,  # type: bool
        req_tracker,  # type: RequirementTracker
        session,  # type: PipSession
        progress_bar,  # type: str
        finder,  # type: PackageFinder
        require_hashes,  # type: bool
        use_user_site,  # type: bool
        lazy_wheel,  # type: bool
    ):
        # type: (...) -> None
        super(RequirementPreparer, self).__init__()

        self.src_dir = src_dir
        self.build_dir = build_dir
        self.req_tracker = req_tracker
        self._session = session
        self._download = Downloader(session, progress_bar)
        self._batch_download = BatchDownloader(session, progress_bar)
        self.finder = finder

        # Where still-packed archives should be written to. If None, they are
        # not saved, and are deleted immediately after unpacking.
        self.download_dir = download_dir

        # Is build isolation allowed?
        self.build_isolation = build_isolation

        # Should hash-checking be required?
        self.require_hashes = require_hashes

        # Should install in user site-packages?
        self.use_user_site = use_user_site

        # Should wheels be downloaded lazily?
        self.use_lazy_wheel = lazy_wheel

        # Memoized downloaded files, as mapping of url: (path, mime type)
        self._downloaded = {}  # type: Dict[str, Tuple[str, str]]

        # Previous "header" printed for a link-based InstallRequirement
        self._previous_requirement_header = ("", "")

    def _log_preparing_link(self, req):
        # type: (InstallRequirement) -> None
        """Provide context for the requirement being prepared."""
        if req.link.is_file and not req.original_link_is_in_wheel_cache:
            message = "Processing %s"
            information = str(display_path(req.link.file_path))
        else:
            message = "Collecting %s"
            information = str(req.req or req)

        if (message, information) != self._previous_requirement_header:
            self._previous_requirement_header = (message, information)
            logger.info(message, information)

        if req.original_link_is_in_wheel_cache:
            with indent_log():
                logger.info("Using cached %s", req.link.filename)

    def _ensure_link_req_src_dir(self, req, parallel_builds):
        # type: (InstallRequirement, bool) -> None
        """Ensure source_dir of a linked InstallRequirement."""
        # Since source_dir is only set for editable requirements.
        if req.link.is_wheel:
            # We don't need to unpack wheels, so no need for a source
            # directory.
            return
        assert req.source_dir is None
        # We always delete unpacked sdists after pip runs.
        req.ensure_has_source_dir(
            self.build_dir,
            autodelete=True,
            parallel_builds=parallel_builds,
        )

        # If a checkout exists, it's unwise to keep going.  version
        # inconsistencies are logged later, but do not fail the
        # installation.
        # FIXME: this won't upgrade when there's an existing
        # package unpacked in `req.source_dir`
        if os.path.exists(os.path.join(req.source_dir, 'setup.py')):
            raise PreviousBuildDirError(
                "pip can't proceed with requirements '{}' due to a"
                "pre-existing build directory ({}). This is likely "
                "due to a previous installation that failed . pip is "
                "being responsible and not assuming it can delete this. "
                "Please delete it and try again.".format(req, req.source_dir)
            )

    def _get_linked_req_hashes(self, req):
        # type: (InstallRequirement) -> Hashes
        # By the time this is called, the requirement's link should have
        # been checked so we can tell what kind of requirements req is
        # and raise some more informative errors than otherwise.
        # (For example, we can raise VcsHashUnsupported for a VCS URL
        # rather than HashMissing.)
        if not self.require_hashes:
            return req.hashes(trust_internet=True)

        # We could check these first 2 conditions inside unpack_url
        # and save repetition of conditions, but then we would
        # report less-useful error messages for unhashable
        # requirements, complaining that there's no hash provided.
        if req.link.is_vcs:
            raise VcsHashUnsupported()
        if req.link.is_existing_dir():
            raise DirectoryUrlHashUnsupported()

        # Unpinned packages are asking for trouble when a new version
        # is uploaded.  This isn't a security check, but it saves users
        # a surprising hash mismatch in the future.
        # file:/// URLs aren't pinnable, so don't complain about them
        # not being pinned.
        if req.original_link is None and not req.is_pinned:
            raise HashUnpinned()

        # If known-good hashes are missing for this requirement,
        # shim it with a facade object that will provoke hash
        # computation and then raise a HashMissing exception
        # showing the user what the hash should be.
        return req.hashes(trust_internet=False) or MissingHashes()

    def _fetch_metadata_using_lazy_wheel(self, link):
        # type: (Link) -> Optional[Distribution]
        """Fetch metadata using lazy wheel, if possible."""
        if not self.use_lazy_wheel:
            return None
        if self.require_hashes:
            logger.debug('Lazy wheel is not used as hash checking is required')
            return None
        if link.is_file or not link.is_wheel:
            logger.debug(
                'Lazy wheel is not used as '
                '%r does not points to a remote wheel',
                link,
            )
            return None

        wheel = Wheel(link.filename)
        name = canonicalize_name(wheel.name)
        logger.info(
            'Obtaining dependency information from %s %s',
            name, wheel.version,
        )
        url = link.url.split('#', 1)[0]
        try:
            return dist_from_wheel_url(name, url, self._session)
        except HTTPRangeRequestUnsupported:
            logger.debug('%s does not support range requests', url)
            return None

    def prepare_linked_requirement(self, req, parallel_builds=False):
        # type: (InstallRequirement, bool) -> Distribution
        """Prepare a requirement to be obtained from req.link."""
        assert req.link
        link = req.link
        self._log_preparing_link(req)
        with indent_log():
            # Check if the relevant file is already available
            # in the download directory
            file_path = None
            if self.download_dir is not None and link.is_wheel:
                hashes = self._get_linked_req_hashes(req)
                file_path = _check_download_dir(req.link, self.download_dir, hashes)

            if file_path is not None:
                # The file is already available, so mark it as downloaded
                self._downloaded[req.link.url] = file_path, None
            else:
                # The file is not available, attempt to fetch only metadata
                wheel_dist = self._fetch_metadata_using_lazy_wheel(link)
                if wheel_dist is not None:
                    req.needs_more_preparation = True
                    return wheel_dist

            # None of the optimizations worked, fully prepare the requirement
            return self._prepare_linked_requirement(req, parallel_builds)

    def prepare_linked_requirements_more(self, reqs, parallel_builds=False):
        # type: (Iterable[InstallRequirement], bool) -> None
        """Prepare a linked requirement more, if needed."""
        reqs = [req for req in reqs if req.needs_more_preparation]
        links = [req.link for req in reqs]

        # Let's download to a temporary directory.
        tmpdir = TempDirectory(kind="unpack", globally_managed=True).path
        self._downloaded.update(self._batch_download(links, tmpdir))
        for req in reqs:
            self._prepare_linked_requirement(req, parallel_builds)

    def _prepare_linked_requirement(self, req, parallel_builds):
        # type: (InstallRequirement, bool) -> Distribution
        assert req.link
        link = req.link

        self._ensure_link_req_src_dir(req, parallel_builds)
        hashes = self._get_linked_req_hashes(req)
        if link.url not in self._downloaded:
            try:
                local_file = unpack_url(
                    link, req.source_dir, self._download,
                    self.download_dir, hashes,
                )
            except NetworkConnectionError as exc:
                raise InstallationError(
                    'Could not install requirement {} because of HTTP '
                    'error {} for URL {}'.format(req, exc, link)
                )
        else:
            file_path, content_type = self._downloaded[link.url]
            if hashes:
                hashes.check_against_path(file_path)
            local_file = File(file_path, content_type)

        # For use in later processing,
        # preserve the file path on the requirement.
        if local_file:
            req.local_file_path = local_file.path

        dist = _get_prepared_distribution(
            req, self.req_tracker, self.finder, self.build_isolation,
        )
        return dist

    def save_linked_requirement(self, req):
        # type: (InstallRequirement) -> None
        assert self.download_dir is not None
        assert req.link is not None
        link = req.link
        if link.is_vcs or (link.is_existing_dir() and req.editable):
            # Make a .zip of the source_dir we already created.
            req.archive(self.download_dir)
            return

        if link.is_existing_dir():
            logger.debug(
                'Not copying link to destination directory '
                'since it is a directory: %s', link,
            )
            return
        if req.local_file_path is None:
            # No distribution was downloaded for this requirement.
            return

        download_location = os.path.join(self.download_dir, link.filename)
        if not os.path.exists(download_location):
            shutil.copy(req.local_file_path, download_location)
            download_path = display_path(download_location)
            logger.info('Saved %s', download_path)

    def prepare_editable_requirement(
        self,
        req,  # type: InstallRequirement
    ):
        # type: (...) -> Distribution
        """Prepare an editable requirement
        """
        assert req.editable, "cannot prepare a non-editable req as editable"

        logger.info('Obtaining %s', req)

        with indent_log():
            if self.require_hashes:
                raise InstallationError(
                    'The editable requirement {} cannot be installed when '
                    'requiring hashes, because there is no single file to '
                    'hash.'.format(req)
                )
            req.ensure_has_source_dir(self.src_dir)
            req.update_editable(self.download_dir is None)

            dist = _get_prepared_distribution(
                req, self.req_tracker, self.finder, self.build_isolation,
            )

            req.check_if_exists(self.use_user_site)

        return dist

    def prepare_installed_requirement(
        self,
        req,  # type: InstallRequirement
        skip_reason  # type: str
    ):
        # type: (...) -> Distribution
        """Prepare an already-installed requirement
        """
        assert req.satisfied_by, "req should have been satisfied but isn't"
        assert skip_reason is not None, (
            "did not get skip reason skipped but req.satisfied_by "
            "is set to {}".format(req.satisfied_by)
        )
        logger.info(
            'Requirement %s: %s (%s)',
            skip_reason, req, req.satisfied_by.version
        )
        with indent_log():
            if self.require_hashes:
                logger.debug(
                    'Since it is already installed, we are trusting this '
                    'package without checking its hash. To ensure a '
                    'completely repeatable environment, install into an '
                    'empty virtualenv.'
                )
            return InstalledDistribution(req).get_pkg_resources_distribution()
