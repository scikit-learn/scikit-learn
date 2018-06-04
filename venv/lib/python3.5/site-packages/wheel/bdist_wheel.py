"""
Create a wheel (.whl) distribution.

A wheel is a built archive format.
"""

import csv
import hashlib
import os
import subprocess
import shutil
import sys
import re
from email.generator import Generator
from distutils.core import Command
from distutils.sysconfig import get_python_version
from distutils import log as logger
from shutil import rmtree

import pkg_resources

from .pep425tags import get_abbr_impl, get_impl_ver, get_abi_tag, get_platform
from .util import native, open_for_csv
from .archive import archive_wheelfile
from .pkginfo import write_pkg_info
from .metadata import pkginfo_to_metadata
from . import pep425tags
from . import __version__ as wheel_version


safe_name = pkg_resources.safe_name
safe_version = pkg_resources.safe_version

PY_LIMITED_API_PATTERN = r'cp3\d'


def safer_name(name):
    return safe_name(name).replace('-', '_')


def safer_version(version):
    return safe_version(version).replace('-', '_')


class bdist_wheel(Command):

    description = 'create a wheel distribution'

    user_options = [('bdist-dir=', 'b',
                     "temporary directory for creating the distribution"),
                    ('plat-name=', 'p',
                     "platform name to embed in generated filenames "
                     "(default: %s)" % get_platform()),
                    ('keep-temp', 'k',
                     "keep the pseudo-installation tree around after " +
                     "creating the distribution archive"),
                    ('dist-dir=', 'd',
                     "directory to put final built distributions in"),
                    ('skip-build', None,
                     "skip rebuilding everything (for testing/debugging)"),
                    ('relative', None,
                     "build the archive using relative paths"
                     "(default: false)"),
                    ('owner=', 'u',
                     "Owner name used when creating a tar file"
                     " [default: current user]"),
                    ('group=', 'g',
                     "Group name used when creating a tar file"
                     " [default: current group]"),
                    ('universal', None,
                     "make a universal wheel"
                     " (default: false)"),
                    ('python-tag=', None,
                     "Python implementation compatibility tag"
                     " (default: py%s)" % get_impl_ver()[0]),
                    ('build-number=', None,
                     "Build number for this particular version. "
                     "As specified in PEP-0427, this must start with a digit. "
                     "[default: None]"),
                    ('py-limited-api=', None,
                     "Python tag (cp32|cp33|cpNN) for abi3 wheel tag"
                     " (default: false)"),
                    ]

    boolean_options = ['keep-temp', 'skip-build', 'relative', 'universal']

    def initialize_options(self):
        self.bdist_dir = None
        self.data_dir = None
        self.plat_name = None
        self.plat_tag = None
        self.format = 'zip'
        self.keep_temp = False
        self.dist_dir = None
        self.distinfo_dir = None
        self.egginfo_dir = None
        self.root_is_pure = None
        self.skip_build = None
        self.relative = False
        self.owner = None
        self.group = None
        self.universal = False
        self.python_tag = 'py' + get_impl_ver()[0]
        self.build_number = None
        self.py_limited_api = False
        self.plat_name_supplied = False

    def finalize_options(self):
        if self.bdist_dir is None:
            bdist_base = self.get_finalized_command('bdist').bdist_base
            self.bdist_dir = os.path.join(bdist_base, 'wheel')

        self.data_dir = self.wheel_dist_name + '.data'
        self.plat_name_supplied = self.plat_name is not None

        need_options = ('dist_dir', 'plat_name', 'skip_build')

        self.set_undefined_options('bdist',
                                   *zip(need_options, need_options))

        self.root_is_pure = not (self.distribution.has_ext_modules()
                                 or self.distribution.has_c_libraries())

        if self.py_limited_api and not re.match(PY_LIMITED_API_PATTERN, self.py_limited_api):
            raise ValueError("py-limited-api must match '%s'" % PY_LIMITED_API_PATTERN)

        # Support legacy [wheel] section for setting universal
        wheel = self.distribution.get_option_dict('wheel')
        if 'universal' in wheel:
            # please don't define this in your global configs
            val = wheel['universal'][1].strip()
            if val.lower() in ('1', 'true', 'yes'):
                self.universal = True

        if self.build_number is not None and not self.build_number[:1].isdigit():
            raise ValueError("Build tag (build-number) must start with a digit.")

    @property
    def wheel_dist_name(self):
        """Return distribution full name with - replaced with _"""
        components = (safer_name(self.distribution.get_name()),
                      safer_version(self.distribution.get_version()))
        if self.build_number:
            components += (self.build_number,)
        return '-'.join(components)

    def get_tag(self):
        # bdist sets self.plat_name if unset, we should only use it for purepy
        # wheels if the user supplied it.
        if self.plat_name_supplied:
            plat_name = self.plat_name
        elif self.root_is_pure:
            plat_name = 'any'
        else:
            plat_name = self.plat_name or get_platform()
            if plat_name in ('linux-x86_64', 'linux_x86_64') and sys.maxsize == 2147483647:
                plat_name = 'linux_i686'
        plat_name = plat_name.replace('-', '_').replace('.', '_')

        if self.root_is_pure:
            if self.universal:
                impl = 'py2.py3'
            else:
                impl = self.python_tag
            tag = (impl, 'none', plat_name)
        else:
            impl_name = get_abbr_impl()
            impl_ver = get_impl_ver()
            impl = impl_name + impl_ver
            # We don't work on CPython 3.1, 3.0.
            if self.py_limited_api and (impl_name + impl_ver).startswith('cp3'):
                impl = self.py_limited_api
                abi_tag = 'abi3'
            else:
                abi_tag = str(get_abi_tag()).lower()
            tag = (impl, abi_tag, plat_name)
            supported_tags = pep425tags.get_supported(
                supplied_platform=plat_name if self.plat_name_supplied else None)
            # XXX switch to this alternate implementation for non-pure:
            if not self.py_limited_api:
                assert tag == supported_tags[0], "%s != %s" % (tag, supported_tags[0])
            assert tag in supported_tags, "would build wheel with unsupported tag {}".format(tag)
        return tag

    def get_archive_basename(self):
        """Return archive name without extension"""

        impl_tag, abi_tag, plat_tag = self.get_tag()

        archive_basename = "%s-%s-%s-%s" % (
            self.wheel_dist_name,
            impl_tag,
            abi_tag,
            plat_tag)
        return archive_basename

    def run(self):
        build_scripts = self.reinitialize_command('build_scripts')
        build_scripts.executable = 'python'

        if not self.skip_build:
            self.run_command('build')

        install = self.reinitialize_command('install',
                                            reinit_subcommands=True)
        install.root = self.bdist_dir
        install.compile = False
        install.skip_build = self.skip_build
        install.warn_dir = False

        # A wheel without setuptools scripts is more cross-platform.
        # Use the (undocumented) `no_ep` option to setuptools'
        # install_scripts command to avoid creating entry point scripts.
        install_scripts = self.reinitialize_command('install_scripts')
        install_scripts.no_ep = True

        # Use a custom scheme for the archive, because we have to decide
        # at installation time which scheme to use.
        for key in ('headers', 'scripts', 'data', 'purelib', 'platlib'):
            setattr(install,
                    'install_' + key,
                    os.path.join(self.data_dir, key))

        basedir_observed = ''

        if os.name == 'nt':
            # win32 barfs if any of these are ''; could be '.'?
            # (distutils.command.install:change_roots bug)
            basedir_observed = os.path.normpath(os.path.join(self.data_dir, '..'))
            self.install_libbase = self.install_lib = basedir_observed

        setattr(install,
                'install_purelib' if self.root_is_pure else 'install_platlib',
                basedir_observed)

        logger.info("installing to %s", self.bdist_dir)

        self.run_command('install')

        archive_basename = self.get_archive_basename()

        pseudoinstall_root = os.path.join(self.dist_dir, archive_basename)
        if not self.relative:
            archive_root = self.bdist_dir
        else:
            archive_root = os.path.join(
                self.bdist_dir,
                self._ensure_relative(install.install_base))

        self.set_undefined_options(
            'install_egg_info', ('target', 'egginfo_dir'))
        self.distinfo_dir = os.path.join(self.bdist_dir,
                                         '%s.dist-info' % self.wheel_dist_name)
        self.egg2dist(self.egginfo_dir,
                      self.distinfo_dir)

        self.write_wheelfile(self.distinfo_dir)

        self.write_record(self.bdist_dir, self.distinfo_dir)

        # Make the archive
        if not os.path.exists(self.dist_dir):
            os.makedirs(self.dist_dir)
        wheel_name = archive_wheelfile(pseudoinstall_root, archive_root)

        # Sign the archive
        if 'WHEEL_TOOL' in os.environ:
            subprocess.call([os.environ['WHEEL_TOOL'], 'sign', wheel_name])

        # Add to 'Distribution.dist_files' so that the "upload" command works
        getattr(self.distribution, 'dist_files', []).append(
            ('bdist_wheel', get_python_version(), wheel_name))

        if not self.keep_temp:
            logger.info('removing %s', self.bdist_dir)
            if not self.dry_run:
                rmtree(self.bdist_dir)

    def write_wheelfile(self, wheelfile_base, generator='bdist_wheel (' + wheel_version + ')'):
        from email.message import Message
        msg = Message()
        msg['Wheel-Version'] = '1.0'  # of the spec
        msg['Generator'] = generator
        msg['Root-Is-Purelib'] = str(self.root_is_pure).lower()
        if self.build_number is not None:
            msg['Build'] = self.build_number

        # Doesn't work for bdist_wininst
        impl_tag, abi_tag, plat_tag = self.get_tag()
        for impl in impl_tag.split('.'):
            for abi in abi_tag.split('.'):
                for plat in plat_tag.split('.'):
                    msg['Tag'] = '-'.join((impl, abi, plat))

        wheelfile_path = os.path.join(wheelfile_base, 'WHEEL')
        logger.info('creating %s', wheelfile_path)
        with open(wheelfile_path, 'w') as f:
            Generator(f, maxheaderlen=0).flatten(msg)

    def _ensure_relative(self, path):
        # copied from dir_util, deleted
        drive, path = os.path.splitdrive(path)
        if path[0:1] == os.sep:
            path = drive + path[1:]
        return path

    def license_file(self):
        """Return license filename from a license-file key in setup.cfg, or None."""
        metadata = self.distribution.get_option_dict('metadata')
        if 'license_file' not in metadata:
            return None
        return metadata['license_file'][1]

    def egg2dist(self, egginfo_path, distinfo_path):
        """Convert an .egg-info directory into a .dist-info directory"""
        def adios(p):
            """Appropriately delete directory, file or link."""
            if os.path.exists(p) and not os.path.islink(p) and os.path.isdir(p):
                shutil.rmtree(p)
            elif os.path.exists(p):
                os.unlink(p)

        adios(distinfo_path)

        if not os.path.exists(egginfo_path):
            # There is no egg-info. This is probably because the egg-info
            # file/directory is not named matching the distribution name used
            # to name the archive file. Check for this case and report
            # accordingly.
            import glob
            pat = os.path.join(os.path.dirname(egginfo_path), '*.egg-info')
            possible = glob.glob(pat)
            err = "Egg metadata expected at %s but not found" % (egginfo_path,)
            if possible:
                alt = os.path.basename(possible[0])
                err += " (%s found - possible misnamed archive file?)" % (alt,)

            raise ValueError(err)

        if os.path.isfile(egginfo_path):
            # .egg-info is a single file
            pkginfo_path = egginfo_path
            pkg_info = pkginfo_to_metadata(egginfo_path, egginfo_path)
            os.mkdir(distinfo_path)
        else:
            # .egg-info is a directory
            pkginfo_path = os.path.join(egginfo_path, 'PKG-INFO')
            pkg_info = pkginfo_to_metadata(egginfo_path, pkginfo_path)

            # ignore common egg metadata that is useless to wheel
            shutil.copytree(egginfo_path, distinfo_path,
                            ignore=lambda x, y: {'PKG-INFO', 'requires.txt', 'SOURCES.txt',
                                                 'not-zip-safe'}
                            )

            # delete dependency_links if it is only whitespace
            dependency_links_path = os.path.join(distinfo_path, 'dependency_links.txt')
            with open(dependency_links_path, 'r') as dependency_links_file:
                dependency_links = dependency_links_file.read().strip()
            if not dependency_links:
                adios(dependency_links_path)

        write_pkg_info(os.path.join(distinfo_path, 'METADATA'), pkg_info)

        # XXX heuristically copy any LICENSE/LICENSE.txt?
        license = self.license_file()
        if license:
            license_filename = 'LICENSE.txt'
            shutil.copy(license, os.path.join(distinfo_path, license_filename))

        adios(egginfo_path)

    def write_record(self, bdist_dir, distinfo_dir):
        from .util import urlsafe_b64encode

        record_path = os.path.join(distinfo_dir, 'RECORD')
        record_relpath = os.path.relpath(record_path, bdist_dir)

        def walk():
            for dir, dirs, files in os.walk(bdist_dir):
                dirs.sort()
                for f in sorted(files):
                    yield os.path.join(dir, f)

        def skip(path):
            """Wheel hashes every possible file."""
            return (path == record_relpath)

        with open_for_csv(record_path, 'w+') as record_file:
            writer = csv.writer(record_file)
            for path in walk():
                relpath = os.path.relpath(path, bdist_dir)
                if skip(relpath):
                    hash = ''
                    size = ''
                else:
                    with open(path, 'rb') as f:
                        data = f.read()
                    digest = hashlib.sha256(data).digest()
                    hash = 'sha256=' + native(urlsafe_b64encode(digest))
                    size = len(data)

                record_path = os.path.relpath(path, bdist_dir).replace(os.path.sep, '/')

                # On Python 2, re-encode the path as UTF-8 from the default file system encoding
                if isinstance(record_path, bytes):
                    record_path = record_path.decode(sys.getfilesystemencoding()).encode('utf-8')

                writer.writerow((record_path, hash, size))
