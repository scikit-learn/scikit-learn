#!python
"""Bootstrap distribute installation

If you want to use setuptools in your package's setup.py, just include this
file in the same directory with it, and add this to the top of your setup.py::

    from distribute_setup import use_setuptools
    use_setuptools()

If you want to require a specific version of setuptools, set a download
mirror, or use an alternate download directory, you can do so by supplying
the appropriate options to ``use_setuptools()``.

This file can also be run as a script to install or upgrade setuptools.
"""
import os
import sys
import time
import fnmatch
import tempfile
import tarfile
from distutils import log

try:
    from site import USER_SITE
except ImportError:
    USER_SITE = None

try:
    import subprocess

    def _python_cmd(*args):
        args = (sys.executable,) + args
        return subprocess.call(args) == 0

except ImportError:
    # will be used for python 2.3
    def _python_cmd(*args):
        args = (sys.executable,) + args
        # quoting arguments if windows
        if sys.platform == 'win32':
            def quote(arg):
                if ' ' in arg:
                    return '"%s"' % arg
                return arg
            args = [quote(arg) for arg in args]
        return os.spawnl(os.P_WAIT, sys.executable, *args) == 0

DEFAULT_VERSION = "0.6.14"
DEFAULT_URL = "http://pypi.python.org/packages/source/d/distribute/"
SETUPTOOLS_FAKED_VERSION = "0.6c11"

SETUPTOOLS_PKG_INFO = """\
Metadata-Version: 1.0
Name: setuptools
Version: %s
Summary: xxxx
Home-page: xxx
Author: xxx
Author-email: xxx
License: xxx
Description: xxx
""" % SETUPTOOLS_FAKED_VERSION


def _install(tarball):
    # extracting the tarball
    tmpdir = tempfile.mkdtemp()
    log.warn('Extracting in %s', tmpdir)
    old_wd = os.getcwd()
    try:
        os.chdir(tmpdir)
        tar = tarfile.open(tarball)
        _extractall(tar)
        tar.close()

        # going in the directory
        subdir = os.path.join(tmpdir, os.listdir(tmpdir)[0])
        os.chdir(subdir)
        log.warn('Now working in %s', subdir)

        # installing
        log.warn('Installing Distribute')
        if not _python_cmd('setup.py', 'install'):
            log.warn('Something went wrong during the installation.')
            log.warn('See the error message above.')
    finally:
        os.chdir(old_wd)


def _build_egg(egg, tarball, to_dir):
    # extracting the tarball
    tmpdir = tempfile.mkdtemp()
    log.warn('Extracting in %s', tmpdir)
    old_wd = os.getcwd()
    try:
        os.chdir(tmpdir)
        tar = tarfile.open(tarball)
        _extractall(tar)
        tar.close()

        # going in the directory
        subdir = os.path.join(tmpdir, os.listdir(tmpdir)[0])
        os.chdir(subdir)
        log.warn('Now working in %s', subdir)

        # building an egg
        log.warn('Building a Distribute egg in %s', to_dir)
        _python_cmd('setup.py', '-q', 'bdist_egg', '--dist-dir', to_dir)

    finally:
        os.chdir(old_wd)
    # returning the result
    log.warn(egg)
    if not os.path.exists(egg):
        raise IOError('Could not build the egg.')


def _do_download(version, download_base, to_dir, download_delay):
    egg = os.path.join(to_dir, 'distribute-%s-py%d.%d.egg'
                       % (version, sys.version_info[0], sys.version_info[1]))
    if not os.path.exists(egg):
        tarball = download_setuptools(version, download_base,
                                      to_dir, download_delay)
        _build_egg(egg, tarball, to_dir)
    sys.path.insert(0, egg)
    import setuptools
    setuptools.bootstrap_install_from = egg


def use_setuptools(version=DEFAULT_VERSION, download_base=DEFAULT_URL,
                   to_dir=os.curdir, download_delay=15, no_fake=True):
    # making sure we use the absolute path
    to_dir = os.path.abspath(to_dir)
    was_imported = 'pkg_resources' in sys.modules or \
        'setuptools' in sys.modules
    try:
        try:
            import pkg_resources
            if not hasattr(pkg_resources, '_distribute'):
                if not no_fake:
                    _fake_setuptools()
                raise ImportError
        except ImportError:
            return _do_download(version, download_base, to_dir, download_delay)
        try:
            pkg_resources.require("distribute>="+version)
            return
        except pkg_resources.VersionConflict:
            e = sys.exc_info()[1]
            if was_imported:
                sys.stderr.write(
                "The required version of distribute (>=%s) is not available,\n"
                "and can't be installed while this script is running. Please\n"
                "install a more recent version first, using\n"
                "'easy_install -U distribute'."
                "\n\n(Currently using %r)\n" % (version, e.args[0]))
                sys.exit(2)
            else:
                del pkg_resources, sys.modules['pkg_resources']    # reload ok
                return _do_download(version, download_base, to_dir,
                                    download_delay)
        except pkg_resources.DistributionNotFound:
            return _do_download(version, download_base, to_dir,
                                download_delay)
    finally:
        if not no_fake:
            _create_fake_setuptools_pkg_info(to_dir)

def download_setuptools(version=DEFAULT_VERSION, download_base=DEFAULT_URL,
                        to_dir=os.curdir, delay=15):
    """Download distribute from a specified location and return its filename

    `version` should be a valid distribute version number that is available
    as an egg for download under the `download_base` URL (which should end
    with a '/'). `to_dir` is the directory where the egg will be downloaded.
    `delay` is the number of seconds to pause before an actual download
    attempt.
    """
    # making sure we use the absolute path
    to_dir = os.path.abspath(to_dir)
    try:
        from urllib.request import urlopen
    except ImportError:
        from urllib2 import urlopen
    tgz_name = "distribute-%s.tar.gz" % version
    url = download_base + tgz_name
    saveto = os.path.join(to_dir, tgz_name)
    src = dst = None
    if not os.path.exists(saveto):  # Avoid repeated downloads
        try:
            log.warn("Downloading %s", url)
            src = urlopen(url)
            # Read/write all in one block, so we don't create a corrupt file
            # if the download is interrupted.
            data = src.read()
            dst = open(saveto, "wb")
            dst.write(data)
        finally:
            if src:
                src.close()
            if dst:
                dst.close()
    return os.path.realpath(saveto)

def _no_sandbox(function):
    def __no_sandbox(*args, **kw):
        try:
            from setuptools.sandbox import DirectorySandbox
            if not hasattr(DirectorySandbox, '_old'):
                def violation(*args):
                    pass
                DirectorySandbox._old = DirectorySandbox._violation
                DirectorySandbox._violation = violation
                patched = True
            else:
                patched = False
        except ImportError:
            patched = False

        try:
            return function(*args, **kw)
        finally:
            if patched:
                DirectorySandbox._violation = DirectorySandbox._old
                del DirectorySandbox._old

    return __no_sandbox

def _patch_file(path, content):
    """Will backup the file then patch it"""
    existing_content = open(path).read()
    if existing_content == content:
        # already patched
        log.warn('Already patched.')
        return False
    log.warn('Patching...')
    _rename_path(path)
    f = open(path, 'w')
    try:
        f.write(content)
    finally:
        f.close()
    return True

_patch_file = _no_sandbox(_patch_file)

def _same_content(path, content):
    return open(path).read() == content

def _rename_path(path):
    new_name = path + '.OLD.%s' % time.time()
    log.warn('Renaming %s into %s', path, new_name)
    os.rename(path, new_name)
    return new_name

def _remove_flat_installation(placeholder):
    if not os.path.isdir(placeholder):
        log.warn('Unkown installation at %s', placeholder)
        return False
    found = False
    for file in os.listdir(placeholder):
        if fnmatch.fnmatch(file, 'setuptools*.egg-info'):
            found = True
            break
    if not found:
        log.warn('Could not locate setuptools*.egg-info')
        return

    log.warn('Removing elements out of the way...')
    pkg_info = os.path.join(placeholder, file)
    if os.path.isdir(pkg_info):
        patched = _patch_egg_dir(pkg_info)
    else:
        patched = _patch_file(pkg_info, SETUPTOOLS_PKG_INFO)

    if not patched:
        log.warn('%s already patched.', pkg_info)
        return False
    # now let's move the files out of the way
    for element in ('setuptools', 'pkg_resources.py', 'site.py'):
        element = os.path.join(placeholder, element)
        if os.path.exists(element):
            _rename_path(element)
        else:
            log.warn('Could not find the %s element of the '
                     'Setuptools distribution', element)
    return True

_remove_flat_installation = _no_sandbox(_remove_flat_installation)

def _after_install(dist):
    log.warn('After install bootstrap.')
    placeholder = dist.get_command_obj('install').install_purelib
    _create_fake_setuptools_pkg_info(placeholder)

def _create_fake_setuptools_pkg_info(placeholder):
    if not placeholder or not os.path.exists(placeholder):
        log.warn('Could not find the install location')
        return
    pyver = '%s.%s' % (sys.version_info[0], sys.version_info[1])
    setuptools_file = 'setuptools-%s-py%s.egg-info' % \
            (SETUPTOOLS_FAKED_VERSION, pyver)
    pkg_info = os.path.join(placeholder, setuptools_file)
    if os.path.exists(pkg_info):
        log.warn('%s already exists', pkg_info)
        return

    log.warn('Creating %s', pkg_info)
    f = open(pkg_info, 'w')
    try:
        f.write(SETUPTOOLS_PKG_INFO)
    finally:
        f.close()

    pth_file = os.path.join(placeholder, 'setuptools.pth')
    log.warn('Creating %s', pth_file)
    f = open(pth_file, 'w')
    try:
        f.write(os.path.join(os.curdir, setuptools_file))
    finally:
        f.close()

_create_fake_setuptools_pkg_info = _no_sandbox(_create_fake_setuptools_pkg_info)

def _patch_egg_dir(path):
    # let's check if it's already patched
    pkg_info = os.path.join(path, 'EGG-INFO', 'PKG-INFO')
    if os.path.exists(pkg_info):
        if _same_content(pkg_info, SETUPTOOLS_PKG_INFO):
            log.warn('%s already patched.', pkg_info)
            return False
    _rename_path(path)
    os.mkdir(path)
    os.mkdir(os.path.join(path, 'EGG-INFO'))
    pkg_info = os.path.join(path, 'EGG-INFO', 'PKG-INFO')
    f = open(pkg_info, 'w')
    try:
        f.write(SETUPTOOLS_PKG_INFO)
    finally:
        f.close()
    return True

_patch_egg_dir = _no_sandbox(_patch_egg_dir)

def _before_install():
    log.warn('Before install bootstrap.')
    _fake_setuptools()


def _under_prefix(location):
    if 'install' not in sys.argv:
        return True
    args = sys.argv[sys.argv.index('install')+1:]
    for index, arg in enumerate(args):
        for option in ('--root', '--prefix'):
            if arg.startswith('%s=' % option):
                top_dir = arg.split('root=')[-1]
                return location.startswith(top_dir)
            elif arg == option:
                if len(args) > index:
                    top_dir = args[index+1]
                    return location.startswith(top_dir)
        if arg == '--user' and USER_SITE is not None:
            return location.startswith(USER_SITE)
    return True


def _fake_setuptools():
    log.warn('Scanning installed packages')
    try:
        import pkg_resources
    except ImportError:
        # we're cool
        log.warn('Setuptools or Distribute does not seem to be installed.')
        return
    ws = pkg_resources.working_set
    try:
        setuptools_dist = ws.find(pkg_resources.Requirement.parse('setuptools',
                                  replacement=False))
    except TypeError:
        # old distribute API
        setuptools_dist = ws.find(pkg_resources.Requirement.parse('setuptools'))

    if setuptools_dist is None:
        log.warn('No setuptools distribution found')
        return
    # detecting if it was already faked
    setuptools_location = setuptools_dist.location
    log.warn('Setuptools installation detected at %s', setuptools_location)

    # if --root or --preix was provided, and if
    # setuptools is not located in them, we don't patch it
    if not _under_prefix(setuptools_location):
        log.warn('Not patching, --root or --prefix is installing Distribute'
                 ' in another location')
        return

    # let's see if its an egg
    if not setuptools_location.endswith('.egg'):
        log.warn('Non-egg installation')
        res = _remove_flat_installation(setuptools_location)
        if not res:
            return
    else:
        log.warn('Egg installation')
        pkg_info = os.path.join(setuptools_location, 'EGG-INFO', 'PKG-INFO')
        if (os.path.exists(pkg_info) and
            _same_content(pkg_info, SETUPTOOLS_PKG_INFO)):
            log.warn('Already patched.')
            return
        log.warn('Patching...')
        # let's create a fake egg replacing setuptools one
        res = _patch_egg_dir(setuptools_location)
        if not res:
            return
    log.warn('Patched done.')
    _relaunch()


def _relaunch():
    log.warn('Relaunching...')
    # we have to relaunch the process
    # pip marker to avoid a relaunch bug
    if sys.argv[:3] == ['-c', 'install', '--single-version-externally-managed']:
        sys.argv[0] = 'setup.py'
    args = [sys.executable] + sys.argv
    sys.exit(subprocess.call(args))


def _extractall(self, path=".", members=None):
    """Extract all members from the archive to the current working
       directory and set owner, modification time and permissions on
       directories afterwards. `path' specifies a different directory
       to extract to. `members' is optional and must be a subset of the
       list returned by getmembers().
    """
    import copy
    import operator
    from tarfile import ExtractError
    directories = []

    if members is None:
        members = self

    for tarinfo in members:
        if tarinfo.isdir():
            # Extract directories with a safe mode.
            directories.append(tarinfo)
            tarinfo = copy.copy(tarinfo)
            tarinfo.mode = 448 # decimal for oct 0700
        self.extract(tarinfo, path)

    # Reverse sort directories.
    if sys.version_info < (2, 4):
        def sorter(dir1, dir2):
            return cmp(dir1.name, dir2.name)
        directories.sort(sorter)
        directories.reverse()
    else:
        directories.sort(key=operator.attrgetter('name'), reverse=True)

    # Set correct owner, mtime and filemode on directories.
    for tarinfo in directories:
        dirpath = os.path.join(path, tarinfo.name)
        try:
            self.chown(tarinfo, dirpath)
            self.utime(tarinfo, dirpath)
            self.chmod(tarinfo, dirpath)
        except ExtractError:
            e = sys.exc_info()[1]
            if self.errorlevel > 1:
                raise
            else:
                self._dbg(1, "tarfile: %s" % e)


def main(argv, version=DEFAULT_VERSION):
    """Install or upgrade setuptools and EasyInstall"""
    tarball = download_setuptools()
    _install(tarball)


if __name__ == '__main__':
    main(sys.argv[1:])
