from __future__ import absolute_import

import io
import logging
import os
import re
import shutil
import sys
import sysconfig
import traceback
import zipfile
from distutils.util import change_root
from email.parser import FeedParser  # type: ignore

from pip._vendor import pkg_resources, pytoml, six
from pip._vendor.packaging import specifiers
from pip._vendor.packaging.markers import Marker
from pip._vendor.packaging.requirements import InvalidRequirement, Requirement
from pip._vendor.packaging.utils import canonicalize_name
from pip._vendor.packaging.version import Version
from pip._vendor.packaging.version import parse as parse_version
from pip._vendor.pkg_resources import RequirementParseError, parse_requirements

from pip._internal import wheel
from pip._internal.build_env import NoOpBuildEnvironment
from pip._internal.compat import native_str
from pip._internal.download import (
    is_archive_file, is_url, path_to_url, url_to_path,
)
from pip._internal.exceptions import InstallationError
from pip._internal.locations import (
    PIP_DELETE_MARKER_FILENAME, running_under_virtualenv,
)
from pip._internal.req.req_uninstall import UninstallPathSet
from pip._internal.utils.hashes import Hashes
from pip._internal.utils.logging import indent_log
from pip._internal.utils.misc import (
    _make_build_dir, ask_path_exists, backup_dir, call_subprocess,
    display_path, dist_in_site_packages, dist_in_usersite, ensure_dir,
    get_installed_version, is_installable_dir, read_text_file, rmtree,
)
from pip._internal.utils.setuptools_build import SETUPTOOLS_SHIM
from pip._internal.utils.temp_dir import TempDirectory
from pip._internal.utils.ui import open_spinner
from pip._internal.vcs import vcs
from pip._internal.wheel import Wheel, move_wheel_files

logger = logging.getLogger(__name__)

operators = specifiers.Specifier._operators.keys()


def _strip_extras(path):
    m = re.match(r'^(.+)(\[[^\]]+\])$', path)
    extras = None
    if m:
        path_no_extras = m.group(1)
        extras = m.group(2)
    else:
        path_no_extras = path

    return path_no_extras, extras


class InstallRequirement(object):
    """
    Represents something that may be installed later on, may have information
    about where to fetch the relavant requirement and also contains logic for
    installing the said requirement.
    """

    def __init__(self, req, comes_from, source_dir=None, editable=False,
                 link=None, update=True, markers=None,
                 isolated=False, options=None, wheel_cache=None,
                 constraint=False, extras=()):
        assert req is None or isinstance(req, Requirement), req
        self.req = req
        self.comes_from = comes_from
        self.constraint = constraint
        if source_dir is not None:
            self.source_dir = os.path.normpath(os.path.abspath(source_dir))
        else:
            self.source_dir = None
        self.editable = editable

        self._wheel_cache = wheel_cache
        if link is not None:
            self.link = self.original_link = link
        else:
            from pip._internal.index import Link
            self.link = self.original_link = req and req.url and Link(req.url)

        if extras:
            self.extras = extras
        elif req:
            self.extras = {
                pkg_resources.safe_extra(extra) for extra in req.extras
            }
        else:
            self.extras = set()
        if markers is not None:
            self.markers = markers
        else:
            self.markers = req and req.marker
        self._egg_info_path = None
        # This holds the pkg_resources.Distribution object if this requirement
        # is already available:
        self.satisfied_by = None
        # This hold the pkg_resources.Distribution object if this requirement
        # conflicts with another installed distribution:
        self.conflicts_with = None
        # Temporary build location
        self._temp_build_dir = TempDirectory(kind="req-build")
        # Used to store the global directory where the _temp_build_dir should
        # have been created. Cf _correct_build_location method.
        self._ideal_build_dir = None
        # True if the editable should be updated:
        self.update = update
        # Set to True after successful installation
        self.install_succeeded = None
        # UninstallPathSet of uninstalled distribution (for possible rollback)
        self.uninstalled_pathset = None
        self.options = options if options else {}
        # Set to True after successful preparation of this requirement
        self.prepared = False
        self.is_direct = False

        self.isolated = isolated
        self.build_env = NoOpBuildEnvironment()

    # Constructors
    #   TODO: Move these out of this class into custom methods.
    @classmethod
    def from_editable(cls, editable_req, comes_from=None, isolated=False,
                      options=None, wheel_cache=None, constraint=False):
        from pip._internal.index import Link

        name, url, extras_override = parse_editable(editable_req)
        if url.startswith('file:'):
            source_dir = url_to_path(url)
        else:
            source_dir = None

        if name is not None:
            try:
                req = Requirement(name)
            except InvalidRequirement:
                raise InstallationError("Invalid requirement: '%s'" % name)
        else:
            req = None
        return cls(
            req, comes_from, source_dir=source_dir,
            editable=True,
            link=Link(url),
            constraint=constraint,
            isolated=isolated,
            options=options if options else {},
            wheel_cache=wheel_cache,
            extras=extras_override or (),
        )

    @classmethod
    def from_req(cls, req, comes_from=None, isolated=False, wheel_cache=None):
        try:
            req = Requirement(req)
        except InvalidRequirement:
            raise InstallationError("Invalid requirement: '%s'" % req)
        if req.url:
            raise InstallationError(
                "Direct url requirement (like %s) are not allowed for "
                "dependencies" % req
            )
        return cls(req, comes_from, isolated=isolated, wheel_cache=wheel_cache)

    @classmethod
    def from_line(
            cls, name, comes_from=None, isolated=False, options=None,
            wheel_cache=None, constraint=False):
        """Creates an InstallRequirement from a name, which might be a
        requirement, directory containing 'setup.py', filename, or URL.
        """
        from pip._internal.index import Link

        if is_url(name):
            marker_sep = '; '
        else:
            marker_sep = ';'
        if marker_sep in name:
            name, markers = name.split(marker_sep, 1)
            markers = markers.strip()
            if not markers:
                markers = None
            else:
                markers = Marker(markers)
        else:
            markers = None
        name = name.strip()
        req = None
        path = os.path.normpath(os.path.abspath(name))
        link = None
        extras = None

        if is_url(name):
            link = Link(name)
        else:
            p, extras = _strip_extras(path)
            looks_like_dir = os.path.isdir(p) and (
                os.path.sep in name or
                (os.path.altsep is not None and os.path.altsep in name) or
                name.startswith('.')
            )
            if looks_like_dir:
                if not is_installable_dir(p):
                    raise InstallationError(
                        "Directory %r is not installable. File 'setup.py' "
                        "not found." % name
                    )
                link = Link(path_to_url(p))
            elif is_archive_file(p):
                if not os.path.isfile(p):
                    logger.warning(
                        'Requirement %r looks like a filename, but the '
                        'file does not exist',
                        name
                    )
                link = Link(path_to_url(p))

        # it's a local file, dir, or url
        if link:
            # Handle relative file URLs
            if link.scheme == 'file' and re.search(r'\.\./', link.url):
                link = Link(
                    path_to_url(os.path.normpath(os.path.abspath(link.path))))
            # wheel file
            if link.is_wheel:
                wheel = Wheel(link.filename)  # can raise InvalidWheelFilename
                req = "%s==%s" % (wheel.name, wheel.version)
            else:
                # set the req to the egg fragment.  when it's not there, this
                # will become an 'unnamed' requirement
                req = link.egg_fragment

        # a requirement specifier
        else:
            req = name

        if extras:
            extras = Requirement("placeholder" + extras.lower()).extras
        else:
            extras = ()
        if req is not None:
            try:
                req = Requirement(req)
            except InvalidRequirement:
                if os.path.sep in req:
                    add_msg = "It looks like a path."
                    add_msg += deduce_helpful_msg(req)
                elif '=' in req and not any(op in req for op in operators):
                    add_msg = "= is not a valid operator. Did you mean == ?"
                else:
                    add_msg = traceback.format_exc()
                raise InstallationError(
                    "Invalid requirement: '%s'\n%s" % (req, add_msg))
        return cls(
            req, comes_from, link=link, markers=markers,
            isolated=isolated,
            options=options if options else {},
            wheel_cache=wheel_cache,
            constraint=constraint,
            extras=extras,
        )

    def __str__(self):
        if self.req:
            s = str(self.req)
            if self.link:
                s += ' from %s' % self.link.url
        else:
            s = self.link.url if self.link else None
        if self.satisfied_by is not None:
            s += ' in %s' % display_path(self.satisfied_by.location)
        if self.comes_from:
            if isinstance(self.comes_from, six.string_types):
                comes_from = self.comes_from
            else:
                comes_from = self.comes_from.from_path()
            if comes_from:
                s += ' (from %s)' % comes_from
        return s

    def __repr__(self):
        return '<%s object: %s editable=%r>' % (
            self.__class__.__name__, str(self), self.editable)

    def populate_link(self, finder, upgrade, require_hashes):
        """Ensure that if a link can be found for this, that it is found.

        Note that self.link may still be None - if Upgrade is False and the
        requirement is already installed.

        If require_hashes is True, don't use the wheel cache, because cached
        wheels, always built locally, have different hashes than the files
        downloaded from the index server and thus throw false hash mismatches.
        Furthermore, cached wheels at present have undeterministic contents due
        to file modification times.
        """
        if self.link is None:
            self.link = finder.find_requirement(self, upgrade)
        if self._wheel_cache is not None and not require_hashes:
            old_link = self.link
            self.link = self._wheel_cache.get(self.link, self.name)
            if old_link != self.link:
                logger.debug('Using cached wheel link: %s', self.link)

    # Things that are valid for all kinds of requirements?
    @property
    def name(self):
        if self.req is None:
            return None
        return native_str(pkg_resources.safe_name(self.req.name))

    @property
    def specifier(self):
        return self.req.specifier

    @property
    def is_pinned(self):
        """Return whether I am pinned to an exact version.

        For example, some-package==1.2 is pinned; some-package>1.2 is not.
        """
        specifiers = self.specifier
        return (len(specifiers) == 1 and
                next(iter(specifiers)).operator in {'==', '==='})

    @property
    def installed_version(self):
        return get_installed_version(self.name)

    def match_markers(self, extras_requested=None):
        if not extras_requested:
            # Provide an extra to safely evaluate the markers
            # without matching any extra
            extras_requested = ('',)
        if self.markers is not None:
            return any(
                self.markers.evaluate({'extra': extra})
                for extra in extras_requested)
        else:
            return True

    @property
    def has_hash_options(self):
        """Return whether any known-good hashes are specified as options.

        These activate --require-hashes mode; hashes specified as part of a
        URL do not.

        """
        return bool(self.options.get('hashes', {}))

    def hashes(self, trust_internet=True):
        """Return a hash-comparer that considers my option- and URL-based
        hashes to be known-good.

        Hashes in URLs--ones embedded in the requirements file, not ones
        downloaded from an index server--are almost peers with ones from
        flags. They satisfy --require-hashes (whether it was implicitly or
        explicitly activated) but do not activate it. md5 and sha224 are not
        allowed in flags, which should nudge people toward good algos. We
        always OR all hashes together, even ones from URLs.

        :param trust_internet: Whether to trust URL-based (#md5=...) hashes
            downloaded from the internet, as by populate_link()

        """
        good_hashes = self.options.get('hashes', {}).copy()
        link = self.link if trust_internet else self.original_link
        if link and link.hash:
            good_hashes.setdefault(link.hash_name, []).append(link.hash)
        return Hashes(good_hashes)

    def from_path(self):
        """Format a nice indicator to show where this "comes from"
        """
        if self.req is None:
            return None
        s = str(self.req)
        if self.comes_from:
            if isinstance(self.comes_from, six.string_types):
                comes_from = self.comes_from
            else:
                comes_from = self.comes_from.from_path()
            if comes_from:
                s += '->' + comes_from
        return s

    def build_location(self, build_dir):
        assert build_dir is not None
        if self._temp_build_dir.path is not None:
            return self._temp_build_dir.path
        if self.req is None:
            # for requirement via a path to a directory: the name of the
            # package is not available yet so we create a temp directory
            # Once run_egg_info will have run, we'll be able
            # to fix it via _correct_build_location
            # Some systems have /tmp as a symlink which confuses custom
            # builds (such as numpy). Thus, we ensure that the real path
            # is returned.
            self._temp_build_dir.create()
            self._ideal_build_dir = build_dir

            return self._temp_build_dir.path
        if self.editable:
            name = self.name.lower()
        else:
            name = self.name
        # FIXME: Is there a better place to create the build_dir? (hg and bzr
        # need this)
        if not os.path.exists(build_dir):
            logger.debug('Creating directory %s', build_dir)
            _make_build_dir(build_dir)
        return os.path.join(build_dir, name)

    def _correct_build_location(self):
        """Move self._temp_build_dir to self._ideal_build_dir/self.req.name

        For some requirements (e.g. a path to a directory), the name of the
        package is not available until we run egg_info, so the build_location
        will return a temporary directory and store the _ideal_build_dir.

        This is only called by self.egg_info_path to fix the temporary build
        directory.
        """
        if self.source_dir is not None:
            return
        assert self.req is not None
        assert self._temp_build_dir.path
        assert self._ideal_build_dir.path
        old_location = self._temp_build_dir.path
        self._temp_build_dir.path = None

        new_location = self.build_location(self._ideal_build_dir)
        if os.path.exists(new_location):
            raise InstallationError(
                'A package already exists in %s; please remove it to continue'
                % display_path(new_location))
        logger.debug(
            'Moving package %s from %s to new location %s',
            self, display_path(old_location), display_path(new_location),
        )
        shutil.move(old_location, new_location)
        self._temp_build_dir.path = new_location
        self._ideal_build_dir = None
        self.source_dir = os.path.normpath(os.path.abspath(new_location))
        self._egg_info_path = None

    def remove_temporary_source(self):
        """Remove the source files from this requirement, if they are marked
        for deletion"""
        if self.source_dir and os.path.exists(
                os.path.join(self.source_dir, PIP_DELETE_MARKER_FILENAME)):
            logger.debug('Removing source in %s', self.source_dir)
            rmtree(self.source_dir)
        self.source_dir = None
        self._temp_build_dir.cleanup()
        self.build_env.cleanup()

    def check_if_exists(self, use_user_site):
        """Find an installed distribution that satisfies or conflicts
        with this requirement, and set self.satisfied_by or
        self.conflicts_with appropriately.
        """
        if self.req is None:
            return False
        try:
            # get_distribution() will resolve the entire list of requirements
            # anyway, and we've already determined that we need the requirement
            # in question, so strip the marker so that we don't try to
            # evaluate it.
            no_marker = Requirement(str(self.req))
            no_marker.marker = None
            self.satisfied_by = pkg_resources.get_distribution(str(no_marker))
            if self.editable and self.satisfied_by:
                self.conflicts_with = self.satisfied_by
                # when installing editables, nothing pre-existing should ever
                # satisfy
                self.satisfied_by = None
                return True
        except pkg_resources.DistributionNotFound:
            return False
        except pkg_resources.VersionConflict:
            existing_dist = pkg_resources.get_distribution(
                self.req.name
            )
            if use_user_site:
                if dist_in_usersite(existing_dist):
                    self.conflicts_with = existing_dist
                elif (running_under_virtualenv() and
                        dist_in_site_packages(existing_dist)):
                    raise InstallationError(
                        "Will not install to the user site because it will "
                        "lack sys.path precedence to %s in %s" %
                        (existing_dist.project_name, existing_dist.location)
                    )
            else:
                self.conflicts_with = existing_dist
        return True

    # Things valid for wheels
    @property
    def is_wheel(self):
        return self.link and self.link.is_wheel

    def move_wheel_files(self, wheeldir, root=None, home=None, prefix=None,
                         warn_script_location=True, use_user_site=False,
                         pycompile=True):
        move_wheel_files(
            self.name, self.req, wheeldir,
            user=use_user_site,
            home=home,
            root=root,
            prefix=prefix,
            pycompile=pycompile,
            isolated=self.isolated,
            warn_script_location=warn_script_location,
        )

    # Things valid for sdists
    @property
    def setup_py_dir(self):
        return os.path.join(
            self.source_dir,
            self.link and self.link.subdirectory_fragment or '')

    @property
    def setup_py(self):
        assert self.source_dir, "No source dir for %s" % self

        setup_py = os.path.join(self.setup_py_dir, 'setup.py')

        # Python2 __file__ should not be unicode
        if six.PY2 and isinstance(setup_py, six.text_type):
            setup_py = setup_py.encode(sys.getfilesystemencoding())

        return setup_py

    @property
    def pyproject_toml(self):
        assert self.source_dir, "No source dir for %s" % self

        pp_toml = os.path.join(self.setup_py_dir, 'pyproject.toml')

        # Python2 __file__ should not be unicode
        if six.PY2 and isinstance(pp_toml, six.text_type):
            pp_toml = pp_toml.encode(sys.getfilesystemencoding())

        return pp_toml

    def get_pep_518_info(self):
        """Get PEP 518 build-time requirements.

        Returns the list of the packages required to build the project,
        specified as per PEP 518 within the package. If `pyproject.toml` is not
        present, returns None to signify not using the same.
        """
        # If pyproject.toml does not exist, don't do anything.
        if not os.path.isfile(self.pyproject_toml):
            return None

        error_template = (
            "{package} has a pyproject.toml file that does not comply "
            "with PEP 518: {reason}"
        )

        with io.open(self.pyproject_toml, encoding="utf-8") as f:
            pp_toml = pytoml.load(f)

        # If there is no build-system table, just use setuptools and wheel.
        if "build-system" not in pp_toml:
            return ["setuptools", "wheel"]

        # Specifying the build-system table but not the requires key is invalid
        build_system = pp_toml["build-system"]
        if "requires" not in build_system:
            raise InstallationError(
                error_template.format(package=self, reason=(
                    "it has a 'build-system' table but not "
                    "'build-system.requires' which is mandatory in the table"
                ))
            )

        # Error out if it's not a list of strings
        requires = build_system["requires"]
        if not _is_list_of_str(requires):
            raise InstallationError(error_template.format(
                package=self,
                reason="'build-system.requires' is not a list of strings.",
            ))

        return requires

    def run_egg_info(self):
        assert self.source_dir
        if self.name:
            logger.debug(
                'Running setup.py (path:%s) egg_info for package %s',
                self.setup_py, self.name,
            )
        else:
            logger.debug(
                'Running setup.py (path:%s) egg_info for package from %s',
                self.setup_py, self.link,
            )

        with indent_log():
            script = SETUPTOOLS_SHIM % self.setup_py
            base_cmd = [sys.executable, '-c', script]
            if self.isolated:
                base_cmd += ["--no-user-cfg"]
            egg_info_cmd = base_cmd + ['egg_info']
            # We can't put the .egg-info files at the root, because then the
            # source code will be mistaken for an installed egg, causing
            # problems
            if self.editable:
                egg_base_option = []
            else:
                egg_info_dir = os.path.join(self.setup_py_dir, 'pip-egg-info')
                ensure_dir(egg_info_dir)
                egg_base_option = ['--egg-base', 'pip-egg-info']
            with self.build_env:
                call_subprocess(
                    egg_info_cmd + egg_base_option,
                    cwd=self.setup_py_dir,
                    show_stdout=False,
                    command_desc='python setup.py egg_info')

        if not self.req:
            if isinstance(parse_version(self.pkg_info()["Version"]), Version):
                op = "=="
            else:
                op = "==="
            self.req = Requirement(
                "".join([
                    self.pkg_info()["Name"],
                    op,
                    self.pkg_info()["Version"],
                ])
            )
            self._correct_build_location()
        else:
            metadata_name = canonicalize_name(self.pkg_info()["Name"])
            if canonicalize_name(self.req.name) != metadata_name:
                logger.warning(
                    'Running setup.py (path:%s) egg_info for package %s '
                    'produced metadata for project name %s. Fix your '
                    '#egg=%s fragments.',
                    self.setup_py, self.name, metadata_name, self.name
                )
                self.req = Requirement(metadata_name)

    def egg_info_data(self, filename):
        if self.satisfied_by is not None:
            if not self.satisfied_by.has_metadata(filename):
                return None
            return self.satisfied_by.get_metadata(filename)
        assert self.source_dir
        filename = self.egg_info_path(filename)
        if not os.path.exists(filename):
            return None
        data = read_text_file(filename)
        return data

    def egg_info_path(self, filename):
        if self._egg_info_path is None:
            if self.editable:
                base = self.source_dir
            else:
                base = os.path.join(self.setup_py_dir, 'pip-egg-info')
            filenames = os.listdir(base)
            if self.editable:
                filenames = []
                for root, dirs, files in os.walk(base):
                    for dir in vcs.dirnames:
                        if dir in dirs:
                            dirs.remove(dir)
                    # Iterate over a copy of ``dirs``, since mutating
                    # a list while iterating over it can cause trouble.
                    # (See https://github.com/pypa/pip/pull/462.)
                    for dir in list(dirs):
                        # Don't search in anything that looks like a virtualenv
                        # environment
                        if (
                                os.path.lexists(
                                    os.path.join(root, dir, 'bin', 'python')
                                ) or
                                os.path.exists(
                                    os.path.join(
                                        root, dir, 'Scripts', 'Python.exe'
                                    )
                                )):
                            dirs.remove(dir)
                        # Also don't search through tests
                        elif dir == 'test' or dir == 'tests':
                            dirs.remove(dir)
                    filenames.extend([os.path.join(root, dir)
                                      for dir in dirs])
                filenames = [f for f in filenames if f.endswith('.egg-info')]

            if not filenames:
                raise InstallationError(
                    "Files/directories (from %s) not found in %s"
                    % (filename, base)
                )
            # if we have more than one match, we pick the toplevel one.  This
            # can easily be the case if there is a dist folder which contains
            # an extracted tarball for testing purposes.
            if len(filenames) > 1:
                filenames.sort(
                    key=lambda x: x.count(os.path.sep) +
                    (os.path.altsep and x.count(os.path.altsep) or 0)
                )
            self._egg_info_path = os.path.join(base, filenames[0])
        return os.path.join(self._egg_info_path, filename)

    def pkg_info(self):
        p = FeedParser()
        data = self.egg_info_data('PKG-INFO')
        if not data:
            logger.warning(
                'No PKG-INFO file found in %s',
                display_path(self.egg_info_path('PKG-INFO')),
            )
        p.feed(data or '')
        return p.close()

    _requirements_section_re = re.compile(r'\[(.*?)\]')

    def get_dist(self):
        """Return a pkg_resources.Distribution built from self.egg_info_path"""
        egg_info = self.egg_info_path('').rstrip(os.path.sep)
        base_dir = os.path.dirname(egg_info)
        metadata = pkg_resources.PathMetadata(base_dir, egg_info)
        dist_name = os.path.splitext(os.path.basename(egg_info))[0]
        return pkg_resources.Distribution(
            os.path.dirname(egg_info),
            project_name=dist_name,
            metadata=metadata,
        )

    def assert_source_matches_version(self):
        assert self.source_dir
        version = self.pkg_info()['version']
        if self.req.specifier and version not in self.req.specifier:
            logger.warning(
                'Requested %s, but installing version %s',
                self,
                version,
            )
        else:
            logger.debug(
                'Source in %s has version %s, which satisfies requirement %s',
                display_path(self.source_dir),
                version,
                self,
            )

    # For both source distributions and editables
    def ensure_has_source_dir(self, parent_dir):
        """Ensure that a source_dir is set.

        This will create a temporary build dir if the name of the requirement
        isn't known yet.

        :param parent_dir: The ideal pip parent_dir for the source_dir.
            Generally src_dir for editables and build_dir for sdists.
        :return: self.source_dir
        """
        if self.source_dir is None:
            self.source_dir = self.build_location(parent_dir)
        return self.source_dir

    # For editable installations
    def install_editable(self, install_options,
                         global_options=(), prefix=None):
        logger.info('Running setup.py develop for %s', self.name)

        if self.isolated:
            global_options = list(global_options) + ["--no-user-cfg"]

        if prefix:
            prefix_param = ['--prefix={}'.format(prefix)]
            install_options = list(install_options) + prefix_param

        with indent_log():
            # FIXME: should we do --install-headers here too?
            with self.build_env:
                call_subprocess(
                    [
                        sys.executable,
                        '-c',
                        SETUPTOOLS_SHIM % self.setup_py
                    ] +
                    list(global_options) +
                    ['develop', '--no-deps'] +
                    list(install_options),

                    cwd=self.setup_py_dir,
                    show_stdout=False,
                )

        self.install_succeeded = True

    def update_editable(self, obtain=True):
        if not self.link:
            logger.debug(
                "Cannot update repository at %s; repository location is "
                "unknown",
                self.source_dir,
            )
            return
        assert self.editable
        assert self.source_dir
        if self.link.scheme == 'file':
            # Static paths don't get updated
            return
        assert '+' in self.link.url, "bad url: %r" % self.link.url
        if not self.update:
            return
        vc_type, url = self.link.url.split('+', 1)
        backend = vcs.get_backend(vc_type)
        if backend:
            vcs_backend = backend(self.link.url)
            if obtain:
                vcs_backend.obtain(self.source_dir)
            else:
                vcs_backend.export(self.source_dir)
        else:
            assert 0, (
                'Unexpected version control type (in %s): %s'
                % (self.link, vc_type))

    # Top-level Actions
    def uninstall(self, auto_confirm=False, verbose=False,
                  use_user_site=False):
        """
        Uninstall the distribution currently satisfying this requirement.

        Prompts before removing or modifying files unless
        ``auto_confirm`` is True.

        Refuses to delete or modify files outside of ``sys.prefix`` -
        thus uninstallation within a virtual environment can only
        modify that virtual environment, even if the virtualenv is
        linked to global site-packages.

        """
        if not self.check_if_exists(use_user_site):
            logger.warning("Skipping %s as it is not installed.", self.name)
            return
        dist = self.satisfied_by or self.conflicts_with

        uninstalled_pathset = UninstallPathSet.from_dist(dist)
        uninstalled_pathset.remove(auto_confirm, verbose)
        return uninstalled_pathset

    def _clean_zip_name(self, name, prefix):  # only used by archive.
        assert name.startswith(prefix + os.path.sep), (
            "name %r doesn't start with prefix %r" % (name, prefix)
        )
        name = name[len(prefix) + 1:]
        name = name.replace(os.path.sep, '/')
        return name

    # TODO: Investigate if this should be kept in InstallRequirement
    #       Seems to be used only when VCS + downloads
    def archive(self, build_dir):
        assert self.source_dir
        create_archive = True
        archive_name = '%s-%s.zip' % (self.name, self.pkg_info()["version"])
        archive_path = os.path.join(build_dir, archive_name)
        if os.path.exists(archive_path):
            response = ask_path_exists(
                'The file %s exists. (i)gnore, (w)ipe, (b)ackup, (a)bort ' %
                display_path(archive_path), ('i', 'w', 'b', 'a'))
            if response == 'i':
                create_archive = False
            elif response == 'w':
                logger.warning('Deleting %s', display_path(archive_path))
                os.remove(archive_path)
            elif response == 'b':
                dest_file = backup_dir(archive_path)
                logger.warning(
                    'Backing up %s to %s',
                    display_path(archive_path),
                    display_path(dest_file),
                )
                shutil.move(archive_path, dest_file)
            elif response == 'a':
                sys.exit(-1)
        if create_archive:
            zip = zipfile.ZipFile(
                archive_path, 'w', zipfile.ZIP_DEFLATED,
                allowZip64=True
            )
            dir = os.path.normcase(os.path.abspath(self.setup_py_dir))
            for dirpath, dirnames, filenames in os.walk(dir):
                if 'pip-egg-info' in dirnames:
                    dirnames.remove('pip-egg-info')
                for dirname in dirnames:
                    dirname = os.path.join(dirpath, dirname)
                    name = self._clean_zip_name(dirname, dir)
                    zipdir = zipfile.ZipInfo(self.name + '/' + name + '/')
                    zipdir.external_attr = 0x1ED << 16  # 0o755
                    zip.writestr(zipdir, '')
                for filename in filenames:
                    if filename == PIP_DELETE_MARKER_FILENAME:
                        continue
                    filename = os.path.join(dirpath, filename)
                    name = self._clean_zip_name(filename, dir)
                    zip.write(filename, self.name + '/' + name)
            zip.close()
            logger.info('Saved %s', display_path(archive_path))

    def install(self, install_options, global_options=None, root=None,
                home=None, prefix=None, warn_script_location=True,
                use_user_site=False, pycompile=True):
        global_options = global_options if global_options is not None else []
        if self.editable:
            self.install_editable(
                install_options, global_options, prefix=prefix,
            )
            return
        if self.is_wheel:
            version = wheel.wheel_version(self.source_dir)
            wheel.check_compatibility(version, self.name)

            self.move_wheel_files(
                self.source_dir, root=root, prefix=prefix, home=home,
                warn_script_location=warn_script_location,
                use_user_site=use_user_site, pycompile=pycompile,
            )
            self.install_succeeded = True
            return

        # Extend the list of global and install options passed on to
        # the setup.py call with the ones from the requirements file.
        # Options specified in requirements file override those
        # specified on the command line, since the last option given
        # to setup.py is the one that is used.
        global_options = list(global_options) + \
            self.options.get('global_options', [])
        install_options = list(install_options) + \
            self.options.get('install_options', [])

        if self.isolated:
            global_options = global_options + ["--no-user-cfg"]

        with TempDirectory(kind="record") as temp_dir:
            record_filename = os.path.join(temp_dir.path, 'install-record.txt')
            install_args = self.get_install_args(
                global_options, record_filename, root, prefix, pycompile,
            )
            msg = 'Running setup.py install for %s' % (self.name,)
            with open_spinner(msg) as spinner:
                with indent_log():
                    with self.build_env:
                        call_subprocess(
                            install_args + install_options,
                            cwd=self.setup_py_dir,
                            show_stdout=False,
                            spinner=spinner,
                        )

            if not os.path.exists(record_filename):
                logger.debug('Record file %s not found', record_filename)
                return
            self.install_succeeded = True

            def prepend_root(path):
                if root is None or not os.path.isabs(path):
                    return path
                else:
                    return change_root(root, path)

            with open(record_filename) as f:
                for line in f:
                    directory = os.path.dirname(line)
                    if directory.endswith('.egg-info'):
                        egg_info_dir = prepend_root(directory)
                        break
                else:
                    logger.warning(
                        'Could not find .egg-info directory in install record'
                        ' for %s',
                        self,
                    )
                    # FIXME: put the record somewhere
                    # FIXME: should this be an error?
                    return
            new_lines = []
            with open(record_filename) as f:
                for line in f:
                    filename = line.strip()
                    if os.path.isdir(filename):
                        filename += os.path.sep
                    new_lines.append(
                        os.path.relpath(prepend_root(filename), egg_info_dir)
                    )
            new_lines.sort()
            ensure_dir(egg_info_dir)
            inst_files_path = os.path.join(egg_info_dir, 'installed-files.txt')
            with open(inst_files_path, 'w') as f:
                f.write('\n'.join(new_lines) + '\n')

    def get_install_args(self, global_options, record_filename, root, prefix,
                         pycompile):
        install_args = [sys.executable, "-u"]
        install_args.append('-c')
        install_args.append(SETUPTOOLS_SHIM % self.setup_py)
        install_args += list(global_options) + \
            ['install', '--record', record_filename]
        install_args += ['--single-version-externally-managed']

        if root is not None:
            install_args += ['--root', root]
        if prefix is not None:
            install_args += ['--prefix', prefix]

        if pycompile:
            install_args += ["--compile"]
        else:
            install_args += ["--no-compile"]

        if running_under_virtualenv():
            py_ver_str = 'python' + sysconfig.get_python_version()
            install_args += ['--install-headers',
                             os.path.join(sys.prefix, 'include', 'site',
                                          py_ver_str, self.name)]

        return install_args


def parse_editable(editable_req):
    """Parses an editable requirement into:
        - a requirement name
        - an URL
        - extras
        - editable options
    Accepted requirements:
        svn+http://blahblah@rev#egg=Foobar[baz]&subdirectory=version_subdir
        .[some_extra]
    """

    from pip._internal.index import Link

    url = editable_req

    # If a file path is specified with extras, strip off the extras.
    url_no_extras, extras = _strip_extras(url)

    if os.path.isdir(url_no_extras):
        if not os.path.exists(os.path.join(url_no_extras, 'setup.py')):
            raise InstallationError(
                "Directory %r is not installable. File 'setup.py' not found." %
                url_no_extras
            )
        # Treating it as code that has already been checked out
        url_no_extras = path_to_url(url_no_extras)

    if url_no_extras.lower().startswith('file:'):
        package_name = Link(url_no_extras).egg_fragment
        if extras:
            return (
                package_name,
                url_no_extras,
                Requirement("placeholder" + extras.lower()).extras,
            )
        else:
            return package_name, url_no_extras, None

    for version_control in vcs:
        if url.lower().startswith('%s:' % version_control):
            url = '%s+%s' % (version_control, url)
            break

    if '+' not in url:
        raise InstallationError(
            '%s should either be a path to a local project or a VCS url '
            'beginning with svn+, git+, hg+, or bzr+' %
            editable_req
        )

    vc_type = url.split('+', 1)[0].lower()

    if not vcs.get_backend(vc_type):
        error_message = 'For --editable=%s only ' % editable_req + \
            ', '.join([backend.name + '+URL' for backend in vcs.backends]) + \
            ' is currently supported'
        raise InstallationError(error_message)

    package_name = Link(url).egg_fragment
    if not package_name:
        raise InstallationError(
            "Could not detect requirement name for '%s', please specify one "
            "with #egg=your_package_name" % editable_req
        )
    return package_name, url, None


def deduce_helpful_msg(req):
    """Returns helpful msg in case requirements file does not exist,
    or cannot be parsed.

    :params req: Requirements file path
    """
    msg = ""
    if os.path.exists(req):
        msg = " It does exist."
        # Try to parse and check if it is a requirements file.
        try:
            with open(req, 'r') as fp:
                # parse first line only
                next(parse_requirements(fp.read()))
                msg += " The argument you provided " + \
                    "(%s) appears to be a" % (req) + \
                    " requirements file. If that is the" + \
                    " case, use the '-r' flag to install" + \
                    " the packages specified within it."
        except RequirementParseError:
            logger.debug("Cannot parse '%s' as requirements \
            file" % (req), exc_info=1)
    else:
        msg += " File '%s' does not exist." % (req)
    return msg


def _is_list_of_str(obj):
    return (
        isinstance(obj, list) and
        all(isinstance(item, six.string_types) for item in obj)
    )
