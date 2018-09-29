"""
Support for installing and building the "wheel" binary package format.
"""
from __future__ import absolute_import

import collections
import compileall
import csv
import hashlib
import logging
import os.path
import re
import shutil
import stat
import sys
import warnings
from base64 import urlsafe_b64encode
from email.parser import Parser

from pip._vendor import pkg_resources
from pip._vendor.distlib.scripts import ScriptMaker
from pip._vendor.packaging.utils import canonicalize_name
from pip._vendor.six import StringIO

from pip._internal import pep425tags
from pip._internal.download import path_to_url, unpack_url
from pip._internal.exceptions import (
    InstallationError, InvalidWheelFilename, UnsupportedWheel,
)
from pip._internal.locations import (
    PIP_DELETE_MARKER_FILENAME, distutils_scheme,
)
from pip._internal.utils.logging import indent_log
from pip._internal.utils.misc import (
    call_subprocess, captured_stdout, ensure_dir, read_chunks,
)
from pip._internal.utils.setuptools_build import SETUPTOOLS_SHIM
from pip._internal.utils.temp_dir import TempDirectory
from pip._internal.utils.typing import MYPY_CHECK_RUNNING
from pip._internal.utils.ui import open_spinner

if MYPY_CHECK_RUNNING:
    from typing import Dict, List, Optional  # noqa: F401

wheel_ext = '.whl'

VERSION_COMPATIBLE = (1, 0)


logger = logging.getLogger(__name__)


def rehash(path, blocksize=1 << 20):
    """Return (hash, length) for path using hashlib.sha256()"""
    h = hashlib.sha256()
    length = 0
    with open(path, 'rb') as f:
        for block in read_chunks(f, size=blocksize):
            length += len(block)
            h.update(block)
    digest = 'sha256=' + urlsafe_b64encode(
        h.digest()
    ).decode('latin1').rstrip('=')
    return (digest, length)


def open_for_csv(name, mode):
    if sys.version_info[0] < 3:
        nl = {}
        bin = 'b'
    else:
        nl = {'newline': ''}
        bin = ''
    return open(name, mode + bin, **nl)


def fix_script(path):
    """Replace #!python with #!/path/to/python
    Return True if file was changed."""
    # XXX RECORD hashes will need to be updated
    if os.path.isfile(path):
        with open(path, 'rb') as script:
            firstline = script.readline()
            if not firstline.startswith(b'#!python'):
                return False
            exename = sys.executable.encode(sys.getfilesystemencoding())
            firstline = b'#!' + exename + os.linesep.encode("ascii")
            rest = script.read()
        with open(path, 'wb') as script:
            script.write(firstline)
            script.write(rest)
        return True


dist_info_re = re.compile(r"""^(?P<namever>(?P<name>.+?)(-(?P<ver>.+?))?)
                                \.dist-info$""", re.VERBOSE)


def root_is_purelib(name, wheeldir):
    """
    Return True if the extracted wheel in wheeldir should go into purelib.
    """
    name_folded = name.replace("-", "_")
    for item in os.listdir(wheeldir):
        match = dist_info_re.match(item)
        if match and match.group('name') == name_folded:
            with open(os.path.join(wheeldir, item, 'WHEEL')) as wheel:
                for line in wheel:
                    line = line.lower().rstrip()
                    if line == "root-is-purelib: true":
                        return True
    return False


def get_entrypoints(filename):
    if not os.path.exists(filename):
        return {}, {}

    # This is done because you can pass a string to entry_points wrappers which
    # means that they may or may not be valid INI files. The attempt here is to
    # strip leading and trailing whitespace in order to make them valid INI
    # files.
    with open(filename) as fp:
        data = StringIO()
        for line in fp:
            data.write(line.strip())
            data.write("\n")
        data.seek(0)

    # get the entry points and then the script names
    entry_points = pkg_resources.EntryPoint.parse_map(data)
    console = entry_points.get('console_scripts', {})
    gui = entry_points.get('gui_scripts', {})

    def _split_ep(s):
        """get the string representation of EntryPoint, remove space and split
        on '='"""
        return str(s).replace(" ", "").split("=")

    # convert the EntryPoint objects into strings with module:function
    console = dict(_split_ep(v) for v in console.values())
    gui = dict(_split_ep(v) for v in gui.values())
    return console, gui


def message_about_scripts_not_on_PATH(scripts):
    # type: (List[str]) -> Optional[str]
    """Determine if any scripts are not on PATH and format a warning.

    Returns a warning message if one or more scripts are not on PATH,
    otherwise None.
    """
    if not scripts:
        return None

    # Group scripts by the path they were installed in
    grouped_by_dir = collections.defaultdict(set)  # type: Dict[str, set]
    for destfile in scripts:
        parent_dir = os.path.dirname(destfile)
        script_name = os.path.basename(destfile)
        grouped_by_dir[parent_dir].add(script_name)

    # We don't want to warn for directories that are on PATH.
    not_warn_dirs = [
        os.path.normcase(i).rstrip(os.sep) for i in
        os.environ.get("PATH", "").split(os.pathsep)
    ]
    # If an executable sits with sys.executable, we don't warn for it.
    #     This covers the case of venv invocations without activating the venv.
    not_warn_dirs.append(os.path.normcase(os.path.dirname(sys.executable)))
    warn_for = {
        parent_dir: scripts for parent_dir, scripts in grouped_by_dir.items()
        if os.path.normcase(parent_dir) not in not_warn_dirs
    }
    if not warn_for:
        return None

    # Format a message
    msg_lines = []
    for parent_dir, scripts in warn_for.items():
        scripts = sorted(scripts)
        if len(scripts) == 1:
            start_text = "script {} is".format(scripts[0])
        else:
            start_text = "scripts {} are".format(
                ", ".join(scripts[:-1]) + " and " + scripts[-1]
            )

        msg_lines.append(
            "The {} installed in '{}' which is not on PATH."
            .format(start_text, parent_dir)
        )

    last_line_fmt = (
        "Consider adding {} to PATH or, if you prefer "
        "to suppress this warning, use --no-warn-script-location."
    )
    if len(msg_lines) == 1:
        msg_lines.append(last_line_fmt.format("this directory"))
    else:
        msg_lines.append(last_line_fmt.format("these directories"))

    # Returns the formatted multiline message
    return "\n".join(msg_lines)


def move_wheel_files(name, req, wheeldir, user=False, home=None, root=None,
                     pycompile=True, scheme=None, isolated=False, prefix=None,
                     warn_script_location=True):
    """Install a wheel"""

    if not scheme:
        scheme = distutils_scheme(
            name, user=user, home=home, root=root, isolated=isolated,
            prefix=prefix,
        )

    if root_is_purelib(name, wheeldir):
        lib_dir = scheme['purelib']
    else:
        lib_dir = scheme['platlib']

    info_dir = []
    data_dirs = []
    source = wheeldir.rstrip(os.path.sep) + os.path.sep

    # Record details of the files moved
    #   installed = files copied from the wheel to the destination
    #   changed = files changed while installing (scripts #! line typically)
    #   generated = files newly generated during the install (script wrappers)
    installed = {}
    changed = set()
    generated = []

    # Compile all of the pyc files that we're going to be installing
    if pycompile:
        with captured_stdout() as stdout:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                compileall.compile_dir(source, force=True, quiet=True)
        logger.debug(stdout.getvalue())

    def normpath(src, p):
        return os.path.relpath(src, p).replace(os.path.sep, '/')

    def record_installed(srcfile, destfile, modified=False):
        """Map archive RECORD paths to installation RECORD paths."""
        oldpath = normpath(srcfile, wheeldir)
        newpath = normpath(destfile, lib_dir)
        installed[oldpath] = newpath
        if modified:
            changed.add(destfile)

    def clobber(source, dest, is_base, fixer=None, filter=None):
        ensure_dir(dest)  # common for the 'include' path

        for dir, subdirs, files in os.walk(source):
            basedir = dir[len(source):].lstrip(os.path.sep)
            destdir = os.path.join(dest, basedir)
            if is_base and basedir.split(os.path.sep, 1)[0].endswith('.data'):
                continue
            for s in subdirs:
                destsubdir = os.path.join(dest, basedir, s)
                if is_base and basedir == '' and destsubdir.endswith('.data'):
                    data_dirs.append(s)
                    continue
                elif (is_base and
                        s.endswith('.dist-info') and
                        canonicalize_name(s).startswith(
                            canonicalize_name(req.name))):
                    assert not info_dir, ('Multiple .dist-info directories: ' +
                                          destsubdir + ', ' +
                                          ', '.join(info_dir))
                    info_dir.append(destsubdir)
            for f in files:
                # Skip unwanted files
                if filter and filter(f):
                    continue
                srcfile = os.path.join(dir, f)
                destfile = os.path.join(dest, basedir, f)
                # directory creation is lazy and after the file filtering above
                # to ensure we don't install empty dirs; empty dirs can't be
                # uninstalled.
                ensure_dir(destdir)

                # copyfile (called below) truncates the destination if it
                # exists and then writes the new contents. This is fine in most
                # cases, but can cause a segfault if pip has loaded a shared
                # object (e.g. from pyopenssl through its vendored urllib3)
                # Since the shared object is mmap'd an attempt to call a
                # symbol in it will then cause a segfault. Unlinking the file
                # allows writing of new contents while allowing the process to
                # continue to use the old copy.
                if os.path.exists(destfile):
                    os.unlink(destfile)

                # We use copyfile (not move, copy, or copy2) to be extra sure
                # that we are not moving directories over (copyfile fails for
                # directories) as well as to ensure that we are not copying
                # over any metadata because we want more control over what
                # metadata we actually copy over.
                shutil.copyfile(srcfile, destfile)

                # Copy over the metadata for the file, currently this only
                # includes the atime and mtime.
                st = os.stat(srcfile)
                if hasattr(os, "utime"):
                    os.utime(destfile, (st.st_atime, st.st_mtime))

                # If our file is executable, then make our destination file
                # executable.
                if os.access(srcfile, os.X_OK):
                    st = os.stat(srcfile)
                    permissions = (
                        st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
                    )
                    os.chmod(destfile, permissions)

                changed = False
                if fixer:
                    changed = fixer(destfile)
                record_installed(srcfile, destfile, changed)

    clobber(source, lib_dir, True)

    assert info_dir, "%s .dist-info directory not found" % req

    # Get the defined entry points
    ep_file = os.path.join(info_dir[0], 'entry_points.txt')
    console, gui = get_entrypoints(ep_file)

    def is_entrypoint_wrapper(name):
        # EP, EP.exe and EP-script.py are scripts generated for
        # entry point EP by setuptools
        if name.lower().endswith('.exe'):
            matchname = name[:-4]
        elif name.lower().endswith('-script.py'):
            matchname = name[:-10]
        elif name.lower().endswith(".pya"):
            matchname = name[:-4]
        else:
            matchname = name
        # Ignore setuptools-generated scripts
        return (matchname in console or matchname in gui)

    for datadir in data_dirs:
        fixer = None
        filter = None
        for subdir in os.listdir(os.path.join(wheeldir, datadir)):
            fixer = None
            if subdir == 'scripts':
                fixer = fix_script
                filter = is_entrypoint_wrapper
            source = os.path.join(wheeldir, datadir, subdir)
            dest = scheme[subdir]
            clobber(source, dest, False, fixer=fixer, filter=filter)

    maker = ScriptMaker(None, scheme['scripts'])

    # Ensure old scripts are overwritten.
    # See https://github.com/pypa/pip/issues/1800
    maker.clobber = True

    # Ensure we don't generate any variants for scripts because this is almost
    # never what somebody wants.
    # See https://bitbucket.org/pypa/distlib/issue/35/
    maker.variants = {''}

    # This is required because otherwise distlib creates scripts that are not
    # executable.
    # See https://bitbucket.org/pypa/distlib/issue/32/
    maker.set_mode = True

    # Simplify the script and fix the fact that the default script swallows
    # every single stack trace.
    # See https://bitbucket.org/pypa/distlib/issue/34/
    # See https://bitbucket.org/pypa/distlib/issue/33/
    def _get_script_text(entry):
        if entry.suffix is None:
            raise InstallationError(
                "Invalid script entry point: %s for req: %s - A callable "
                "suffix is required. Cf https://packaging.python.org/en/"
                "latest/distributing.html#console-scripts for more "
                "information." % (entry, req)
            )
        return maker.script_template % {
            "module": entry.prefix,
            "import_name": entry.suffix.split(".")[0],
            "func": entry.suffix,
        }

    maker._get_script_text = _get_script_text
    maker.script_template = r"""# -*- coding: utf-8 -*-
import re
import sys

from %(module)s import %(import_name)s

if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0])
    sys.exit(%(func)s())
"""

    # Special case pip and setuptools to generate versioned wrappers
    #
    # The issue is that some projects (specifically, pip and setuptools) use
    # code in setup.py to create "versioned" entry points - pip2.7 on Python
    # 2.7, pip3.3 on Python 3.3, etc. But these entry points are baked into
    # the wheel metadata at build time, and so if the wheel is installed with
    # a *different* version of Python the entry points will be wrong. The
    # correct fix for this is to enhance the metadata to be able to describe
    # such versioned entry points, but that won't happen till Metadata 2.0 is
    # available.
    # In the meantime, projects using versioned entry points will either have
    # incorrect versioned entry points, or they will not be able to distribute
    # "universal" wheels (i.e., they will need a wheel per Python version).
    #
    # Because setuptools and pip are bundled with _ensurepip and virtualenv,
    # we need to use universal wheels. So, as a stopgap until Metadata 2.0, we
    # override the versioned entry points in the wheel and generate the
    # correct ones. This code is purely a short-term measure until Metadata 2.0
    # is available.
    #
    # To add the level of hack in this section of code, in order to support
    # ensurepip this code will look for an ``ENSUREPIP_OPTIONS`` environment
    # variable which will control which version scripts get installed.
    #
    # ENSUREPIP_OPTIONS=altinstall
    #   - Only pipX.Y and easy_install-X.Y will be generated and installed
    # ENSUREPIP_OPTIONS=install
    #   - pipX.Y, pipX, easy_install-X.Y will be generated and installed. Note
    #     that this option is technically if ENSUREPIP_OPTIONS is set and is
    #     not altinstall
    # DEFAULT
    #   - The default behavior is to install pip, pipX, pipX.Y, easy_install
    #     and easy_install-X.Y.
    pip_script = console.pop('pip', None)
    if pip_script:
        if "ENSUREPIP_OPTIONS" not in os.environ:
            spec = 'pip = ' + pip_script
            generated.extend(maker.make(spec))

        if os.environ.get("ENSUREPIP_OPTIONS", "") != "altinstall":
            spec = 'pip%s = %s' % (sys.version[:1], pip_script)
            generated.extend(maker.make(spec))

        spec = 'pip%s = %s' % (sys.version[:3], pip_script)
        generated.extend(maker.make(spec))
        # Delete any other versioned pip entry points
        pip_ep = [k for k in console if re.match(r'pip(\d(\.\d)?)?$', k)]
        for k in pip_ep:
            del console[k]
    easy_install_script = console.pop('easy_install', None)
    if easy_install_script:
        if "ENSUREPIP_OPTIONS" not in os.environ:
            spec = 'easy_install = ' + easy_install_script
            generated.extend(maker.make(spec))

        spec = 'easy_install-%s = %s' % (sys.version[:3], easy_install_script)
        generated.extend(maker.make(spec))
        # Delete any other versioned easy_install entry points
        easy_install_ep = [
            k for k in console if re.match(r'easy_install(-\d\.\d)?$', k)
        ]
        for k in easy_install_ep:
            del console[k]

    # Generate the console and GUI entry points specified in the wheel
    if len(console) > 0:
        generated_console_scripts = maker.make_multiple(
            ['%s = %s' % kv for kv in console.items()]
        )
        generated.extend(generated_console_scripts)

        if warn_script_location:
            msg = message_about_scripts_not_on_PATH(generated_console_scripts)
            if msg is not None:
                logger.warn(msg)

    if len(gui) > 0:
        generated.extend(
            maker.make_multiple(
                ['%s = %s' % kv for kv in gui.items()],
                {'gui': True}
            )
        )

    # Record pip as the installer
    installer = os.path.join(info_dir[0], 'INSTALLER')
    temp_installer = os.path.join(info_dir[0], 'INSTALLER.pip')
    with open(temp_installer, 'wb') as installer_file:
        installer_file.write(b'pip\n')
    shutil.move(temp_installer, installer)
    generated.append(installer)

    # Record details of all files installed
    record = os.path.join(info_dir[0], 'RECORD')
    temp_record = os.path.join(info_dir[0], 'RECORD.pip')
    with open_for_csv(record, 'r') as record_in:
        with open_for_csv(temp_record, 'w+') as record_out:
            reader = csv.reader(record_in)
            writer = csv.writer(record_out)
            for row in reader:
                row[0] = installed.pop(row[0], row[0])
                if row[0] in changed:
                    row[1], row[2] = rehash(row[0])
                writer.writerow(row)
            for f in generated:
                digest, length = rehash(f)
                writer.writerow((normpath(f, lib_dir), digest, length))
            for f in installed:
                writer.writerow((installed[f], '', ''))
    shutil.move(temp_record, record)


def wheel_version(source_dir):
    """
    Return the Wheel-Version of an extracted wheel, if possible.

    Otherwise, return False if we couldn't parse / extract it.
    """
    try:
        dist = [d for d in pkg_resources.find_on_path(None, source_dir)][0]

        wheel_data = dist.get_metadata('WHEEL')
        wheel_data = Parser().parsestr(wheel_data)

        version = wheel_data['Wheel-Version'].strip()
        version = tuple(map(int, version.split('.')))
        return version
    except Exception:
        return False


def check_compatibility(version, name):
    """
    Raises errors or warns if called with an incompatible Wheel-Version.

    Pip should refuse to install a Wheel-Version that's a major series
    ahead of what it's compatible with (e.g 2.0 > 1.1); and warn when
    installing a version only minor version ahead (e.g 1.2 > 1.1).

    version: a 2-tuple representing a Wheel-Version (Major, Minor)
    name: name of wheel or package to raise exception about

    :raises UnsupportedWheel: when an incompatible Wheel-Version is given
    """
    if not version:
        raise UnsupportedWheel(
            "%s is in an unsupported or invalid wheel" % name
        )
    if version[0] > VERSION_COMPATIBLE[0]:
        raise UnsupportedWheel(
            "%s's Wheel-Version (%s) is not compatible with this version "
            "of pip" % (name, '.'.join(map(str, version)))
        )
    elif version > VERSION_COMPATIBLE:
        logger.warning(
            'Installing from a newer Wheel-Version (%s)',
            '.'.join(map(str, version)),
        )


class Wheel(object):
    """A wheel file"""

    # TODO: maybe move the install code into this class

    wheel_file_re = re.compile(
        r"""^(?P<namever>(?P<name>.+?)-(?P<ver>.*?))
        ((-(?P<build>\d[^-]*?))?-(?P<pyver>.+?)-(?P<abi>.+?)-(?P<plat>.+?)
        \.whl|\.dist-info)$""",
        re.VERBOSE
    )

    def __init__(self, filename):
        """
        :raises InvalidWheelFilename: when the filename is invalid for a wheel
        """
        wheel_info = self.wheel_file_re.match(filename)
        if not wheel_info:
            raise InvalidWheelFilename(
                "%s is not a valid wheel filename." % filename
            )
        self.filename = filename
        self.name = wheel_info.group('name').replace('_', '-')
        # we'll assume "_" means "-" due to wheel naming scheme
        # (https://github.com/pypa/pip/issues/1150)
        self.version = wheel_info.group('ver').replace('_', '-')
        self.build_tag = wheel_info.group('build')
        self.pyversions = wheel_info.group('pyver').split('.')
        self.abis = wheel_info.group('abi').split('.')
        self.plats = wheel_info.group('plat').split('.')

        # All the tag combinations from this file
        self.file_tags = {
            (x, y, z) for x in self.pyversions
            for y in self.abis for z in self.plats
        }

    def support_index_min(self, tags=None):
        """
        Return the lowest index that one of the wheel's file_tag combinations
        achieves in the supported_tags list e.g. if there are 8 supported tags,
        and one of the file tags is first in the list, then return 0.  Returns
        None is the wheel is not supported.
        """
        if tags is None:  # for mock
            tags = pep425tags.get_supported()
        indexes = [tags.index(c) for c in self.file_tags if c in tags]
        return min(indexes) if indexes else None

    def supported(self, tags=None):
        """Is this wheel supported on this system?"""
        if tags is None:  # for mock
            tags = pep425tags.get_supported()
        return bool(set(tags).intersection(self.file_tags))


class WheelBuilder(object):
    """Build wheels from a RequirementSet."""

    def __init__(self, finder, preparer, wheel_cache,
                 build_options=None, global_options=None, no_clean=False):
        self.finder = finder
        self.preparer = preparer
        self.wheel_cache = wheel_cache

        self._wheel_dir = preparer.wheel_download_dir

        self.build_options = build_options or []
        self.global_options = global_options or []
        self.no_clean = no_clean

    def _build_one(self, req, output_dir, python_tag=None):
        """Build one wheel.

        :return: The filename of the built wheel, or None if the build failed.
        """
        # Install build deps into temporary directory (PEP 518)
        with req.build_env:
            return self._build_one_inside_env(req, output_dir,
                                              python_tag=python_tag)

    def _build_one_inside_env(self, req, output_dir, python_tag=None):
        with TempDirectory(kind="wheel") as temp_dir:
            if self.__build_one(req, temp_dir.path, python_tag=python_tag):
                try:
                    wheel_name = os.listdir(temp_dir.path)[0]
                    wheel_path = os.path.join(output_dir, wheel_name)
                    shutil.move(
                        os.path.join(temp_dir.path, wheel_name), wheel_path
                    )
                    logger.info('Stored in directory: %s', output_dir)
                    return wheel_path
                except Exception:
                    pass
            # Ignore return, we can't do anything else useful.
            self._clean_one(req)
            return None

    def _base_setup_args(self, req):
        # NOTE: Eventually, we'd want to also -S to the flags here, when we're
        # isolating. Currently, it breaks Python in virtualenvs, because it
        # relies on site.py to find parts of the standard library outside the
        # virtualenv.
        return [
            sys.executable, '-u', '-c',
            SETUPTOOLS_SHIM % req.setup_py
        ] + list(self.global_options)

    def __build_one(self, req, tempd, python_tag=None):
        base_args = self._base_setup_args(req)

        spin_message = 'Running setup.py bdist_wheel for %s' % (req.name,)
        with open_spinner(spin_message) as spinner:
            logger.debug('Destination directory: %s', tempd)
            wheel_args = base_args + ['bdist_wheel', '-d', tempd] \
                + self.build_options

            if python_tag is not None:
                wheel_args += ["--python-tag", python_tag]

            try:
                call_subprocess(wheel_args, cwd=req.setup_py_dir,
                                show_stdout=False, spinner=spinner)
                return True
            except Exception:
                spinner.finish("error")
                logger.error('Failed building wheel for %s', req.name)
                return False

    def _clean_one(self, req):
        base_args = self._base_setup_args(req)

        logger.info('Running setup.py clean for %s', req.name)
        clean_args = base_args + ['clean', '--all']
        try:
            call_subprocess(clean_args, cwd=req.source_dir, show_stdout=False)
            return True
        except Exception:
            logger.error('Failed cleaning build dir for %s', req.name)
            return False

    def build(self, requirements, session, autobuilding=False):
        """Build wheels.

        :param unpack: If True, replace the sdist we built from with the
            newly built wheel, in preparation for installation.
        :return: True if all the wheels built correctly.
        """
        from pip._internal import index

        building_is_possible = self._wheel_dir or (
            autobuilding and self.wheel_cache.cache_dir
        )
        assert building_is_possible

        buildset = []
        for req in requirements:
            if req.constraint:
                continue
            if req.is_wheel:
                if not autobuilding:
                    logger.info(
                        'Skipping %s, due to already being wheel.', req.name,
                    )
            elif autobuilding and req.editable:
                pass
            elif autobuilding and not req.source_dir:
                pass
            elif autobuilding and req.link and not req.link.is_artifact:
                # VCS checkout. Build wheel just for this run.
                buildset.append((req, True))
            else:
                ephem_cache = False
                if autobuilding:
                    link = req.link
                    base, ext = link.splitext()
                    if index.egg_info_matches(base, None, link) is None:
                        # E.g. local directory. Build wheel just for this run.
                        ephem_cache = True
                    if "binary" not in index.fmt_ctl_formats(
                            self.finder.format_control,
                            canonicalize_name(req.name)):
                        logger.info(
                            "Skipping bdist_wheel for %s, due to binaries "
                            "being disabled for it.", req.name,
                        )
                        continue
                buildset.append((req, ephem_cache))

        if not buildset:
            return True

        # Build the wheels.
        logger.info(
            'Building wheels for collected packages: %s',
            ', '.join([req.name for (req, _) in buildset]),
        )
        _cache = self.wheel_cache  # shorter name
        with indent_log():
            build_success, build_failure = [], []
            for req, ephem in buildset:
                python_tag = None
                if autobuilding:
                    python_tag = pep425tags.implementation_tag
                    if ephem:
                        output_dir = _cache.get_ephem_path_for_link(req.link)
                    else:
                        output_dir = _cache.get_path_for_link(req.link)
                    try:
                        ensure_dir(output_dir)
                    except OSError as e:
                        logger.warning("Building wheel for %s failed: %s",
                                       req.name, e)
                        build_failure.append(req)
                        continue
                else:
                    output_dir = self._wheel_dir
                wheel_file = self._build_one(
                    req, output_dir,
                    python_tag=python_tag,
                )
                if wheel_file:
                    build_success.append(req)
                    if autobuilding:
                        # XXX: This is mildly duplicative with prepare_files,
                        # but not close enough to pull out to a single common
                        # method.
                        # The code below assumes temporary source dirs -
                        # prevent it doing bad things.
                        if req.source_dir and not os.path.exists(os.path.join(
                                req.source_dir, PIP_DELETE_MARKER_FILENAME)):
                            raise AssertionError(
                                "bad source dir - missing marker")
                        # Delete the source we built the wheel from
                        req.remove_temporary_source()
                        # set the build directory again - name is known from
                        # the work prepare_files did.
                        req.source_dir = req.build_location(
                            self.preparer.build_dir
                        )
                        # Update the link for this.
                        req.link = index.Link(path_to_url(wheel_file))
                        assert req.link.is_wheel
                        # extract the wheel into the dir
                        unpack_url(
                            req.link, req.source_dir, None, False,
                            session=session,
                        )
                else:
                    build_failure.append(req)

        # notify success/failure
        if build_success:
            logger.info(
                'Successfully built %s',
                ' '.join([req.name for req in build_success]),
            )
        if build_failure:
            logger.info(
                'Failed to build %s',
                ' '.join([req.name for req in build_failure]),
            )
        # Return True if all builds were successful
        return len(build_failure) == 0
