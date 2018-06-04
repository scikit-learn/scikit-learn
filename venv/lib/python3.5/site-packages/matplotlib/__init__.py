"""
This is an object-oriented plotting library.

A procedural interface is provided by the companion pyplot module,
which may be imported directly, e.g.::

    import matplotlib.pyplot as plt

or using ipython::

    ipython

at your terminal, followed by::

    In [1]: %matplotlib
    In [2]: import matplotlib.pyplot as plt

at the ipython shell prompt.

For the most part, direct use of the object-oriented library is
encouraged when programming; pyplot is primarily for working
interactively.  The
exceptions are the pyplot commands :func:`~matplotlib.pyplot.figure`,
:func:`~matplotlib.pyplot.subplot`,
:func:`~matplotlib.pyplot.subplots`, and
:func:`~pyplot.savefig`, which can greatly simplify scripting.

Modules include:

    :mod:`matplotlib.axes`
        defines the :class:`~matplotlib.axes.Axes` class.  Most pylab
        commands are wrappers for :class:`~matplotlib.axes.Axes`
        methods.  The axes module is the highest level of OO access to
        the library.

    :mod:`matplotlib.figure`
        defines the :class:`~matplotlib.figure.Figure` class.

    :mod:`matplotlib.artist`
        defines the :class:`~matplotlib.artist.Artist` base class for
        all classes that draw things.

    :mod:`matplotlib.lines`
        defines the :class:`~matplotlib.lines.Line2D` class for
        drawing lines and markers

    :mod:`matplotlib.patches`
        defines classes for drawing polygons

    :mod:`matplotlib.text`
        defines the :class:`~matplotlib.text.Text`,
        :class:`~matplotlib.text.TextWithDash`, and
        :class:`~matplotlib.text.Annotate` classes

    :mod:`matplotlib.image`
        defines the :class:`~matplotlib.image.AxesImage` and
        :class:`~matplotlib.image.FigureImage` classes

    :mod:`matplotlib.collections`
        classes for efficient drawing of groups of lines or polygons

    :mod:`matplotlib.colors`
        classes for interpreting color specifications and for making
        colormaps

    :mod:`matplotlib.cm`
        colormaps and the :class:`~matplotlib.image.ScalarMappable`
        mixin class for providing color mapping functionality to other
        classes

    :mod:`matplotlib.ticker`
        classes for calculating tick mark locations and for formatting
        tick labels

    :mod:`matplotlib.backends`
        a subpackage with modules for various gui libraries and output
        formats

The base matplotlib namespace includes:

    :data:`~matplotlib.rcParams`
        a global dictionary of default configuration settings.  It is
        initialized by code which may be overridden by a matplotlibrc
        file.

    :func:`~matplotlib.rc`
        a function for setting groups of rcParams values

    :func:`~matplotlib.use`
        a function for setting the matplotlib backend.  If used, this
        function must be called immediately after importing matplotlib
        for the first time.  In particular, it must be called
        **before** importing pylab (if pylab is imported).

matplotlib was initially written by John D. Hunter (1968-2012) and is now
developed and maintained by a host of others.

Occasionally the internal documentation (python docstrings) will refer
to MATLAB&reg;, a registered trademark of The MathWorks, Inc.

"""
from __future__ import absolute_import, division, print_function

import six

import atexit
from collections import MutableMapping
import contextlib
import distutils.version
import distutils.sysconfig
import functools
import io
import inspect
import itertools
import locale
import logging
import os
import re
import shutil
import stat
import sys
import tempfile
import warnings

# cbook must import matplotlib only within function
# definitions, so it is safe to import from it here.
from . import cbook
from matplotlib.cbook import (
    _backports, mplDeprecation, dedent, get_label, sanitize_sequence)
from matplotlib.compat import subprocess
from matplotlib.rcsetup import defaultParams, validate_backend, cycler

import numpy
from six.moves.urllib.request import urlopen
from six.moves import reload_module as reload

# Get the version from the _version.py versioneer file. For a git checkout,
# this is computed based on the number of commits since the last tag.
from ._version import get_versions
__version__ = str(get_versions()['version'])
del get_versions

_log = logging.getLogger(__name__)

__version__numpy__ = str('1.7.1')  # minimum required numpy version

__bibtex__ = r"""@Article{Hunter:2007,
  Author    = {Hunter, J. D.},
  Title     = {Matplotlib: A 2D graphics environment},
  Journal   = {Computing In Science \& Engineering},
  Volume    = {9},
  Number    = {3},
  Pages     = {90--95},
  abstract  = {Matplotlib is a 2D graphics package used for Python
  for application development, interactive scripting, and
  publication-quality image generation across user
  interfaces and operating systems.},
  publisher = {IEEE COMPUTER SOC},
  year      = 2007
}"""


_python27 = (sys.version_info.major == 2 and sys.version_info.minor >= 7)
_python34 = (sys.version_info.major == 3 and sys.version_info.minor >= 4)
if not (_python27 or _python34):
    raise ImportError("Matplotlib requires Python 2.7 or 3.4 or later")

if _python27:
    _log.addHandler(logging.NullHandler())


def compare_versions(a, b):
    "return True if a is greater than or equal to b"
    if a:
        if six.PY3:
            if isinstance(a, bytes):
                a = a.decode('ascii')
            if isinstance(b, bytes):
                b = b.decode('ascii')
        a = distutils.version.LooseVersion(a)
        b = distutils.version.LooseVersion(b)
        return a >= b
    else:
        return False


try:
    import dateutil
except ImportError:
    raise ImportError("Matplotlib requires dateutil")


if not compare_versions(six.__version__, '1.10'):
    raise ImportError(
        "Matplotlib requires six>=1.10; you have %s" % six.__version__)


try:
    import pyparsing
except ImportError:
    raise ImportError("Matplotlib requires pyparsing")
else:
    if not compare_versions(pyparsing.__version__, '2.0.1'):
        raise ImportError(
            "Matplotlib requires pyparsing>=2.0.1; you have %s"
            % pyparsing.__version__)


if not compare_versions(numpy.__version__, __version__numpy__):
    raise ImportError(
        "Matplotlib requires numpy>=%s; you have %s" % (
            __version__numpy__, numpy.__version__))


if not hasattr(sys, 'argv'):  # for modpython
    sys.argv = [str('modpython')]


def _is_writable_dir(p):
    """
    p is a string pointing to a putative writable dir -- return True p
    is such a string, else False
    """
    return os.access(p, os.W_OK) and os.path.isdir(p)

_verbose_msg = """\
matplotlib.verbose is deprecated;
Command line argument --verbose-LEVEL is deprecated.
This functionality is now provided by the standard
python logging library.  To get more (or less) logging output:
    import logging
    logger = logging.getLogger('matplotlib')
    logger.set_level(logging.INFO)"""


def _set_logger_verbose_level(level_str='silent', file_str='sys.stdout'):
    """
    Use a --verbose-LEVEL level to set the logging level:

    """
    levelmap = {'silent': logging.WARNING, 'helpful': logging.INFO,
                'debug': logging.DEBUG, 'debug-annoying': logging.DEBUG,
                'info': logging.INFO, 'warning': logging.WARNING}
    # Check that current state of logger isn't already more verbose
    # than the requested level.  If it is more verbose, then leave more
    # verbose.
    newlev = levelmap[level_str]
    oldlev = _log.getEffectiveLevel()
    if newlev < oldlev:
        _log.setLevel(newlev)
        std = {
            'sys.stdout': sys.stdout,
            'sys.stderr': sys.stderr,
        }
        if file_str in std:
            fileo = std[file_str]
        else:
            fileo = sys.stdout
            try:
                fileo = open(file_str, 'w')
                # if this fails, we will just write to stdout
            except IOError:
                warnings.warn('could not open log file "{0}"'
                              'for writing.  Check your '
                              'matplotlibrc'.format(file_str))
        console = logging.StreamHandler(fileo)
        console.setLevel(newlev)
        _log.addHandler(console)


def _parse_commandline():
    """
    Check for --verbose-LEVEL type command line arguments and
    set logging level appropriately.
    """

    levels = ('silent', 'helpful', 'debug', 'debug-annoying',
              'info', 'warning')

    for arg in sys.argv[1:]:
        if arg.startswith('--verbose-'):
            level_str = arg[10:]
            # If it doesn't match one of ours, then don't even
            # bother noting it, we are just a 3rd-party library
            # to somebody else's script.
            if level_str in levels:
                _set_logger_verbose_level(level_str)

_parse_commandline()


class Verbose(object):
    """
    A class to handle reporting.  Set the fileo attribute to any file
    instance to handle the output.  Default is sys.stdout
    """
    levels = ('silent', 'helpful', 'debug', 'debug-annoying')
    vald = {level: i for i, level in enumerate(levels)}

    # parse the verbosity from the command line; flags look like
    # --verbose-silent or --verbose-helpful
    _commandLineVerbose = None

    for arg in sys.argv[1:]:
        if not arg.startswith('--verbose-'):
            continue
        level_str = arg[10:]
        # If it doesn't match one of ours, then don't even
        # bother noting it, we are just a 3rd-party library
        # to somebody else's script.
        if level_str in levels:
            _commandLineVerbose = level_str

    @cbook.deprecated("2.2", message=_verbose_msg)
    def __init__(self):
        self.set_level('silent')
        self.fileo = sys.stdout

    @cbook.deprecated("2.2", message=_verbose_msg)
    def set_level(self, level):
        'set the verbosity to one of the Verbose.levels strings'

        if self._commandLineVerbose is not None:
            level = self._commandLineVerbose
        if level not in self.levels:
            warnings.warn('matplotlib: unrecognized --verbose-* string "%s".'
                          ' Legal values are %s' % (level, self.levels))
        else:
            self.level = level

    @cbook.deprecated("2.2", message=_verbose_msg)
    def set_fileo(self, fname):
        std = {
            'sys.stdout': sys.stdout,
            'sys.stderr': sys.stderr,
        }
        if fname in std:
            self.fileo = std[fname]
        else:
            try:
                fileo = open(fname, 'w')
            except IOError:
                raise ValueError('Verbose object could not open log file "{0}"'
                                 ' for writing.\nCheck your matplotlibrc '
                                 'verbose.fileo setting'.format(fname))
            else:
                self.fileo = fileo

    @cbook.deprecated("2.2", message=_verbose_msg)
    def report(self, s, level='helpful'):
        """
        print message s to self.fileo if self.level>=level.  Return
        value indicates whether a message was issued

        """
        if self.ge(level):
            print(s, file=self.fileo)
            return True
        return False

    @cbook.deprecated("2.2", message=_verbose_msg)
    def wrap(self, fmt, func, level='helpful', always=True):
        """
        return a callable function that wraps func and reports it
        output through the verbose handler if current verbosity level
        is higher than level

        if always is True, the report will occur on every function
        call; otherwise only on the first time the function is called
        """
        assert callable(func)

        def wrapper(*args, **kwargs):
            ret = func(*args, **kwargs)

            if (always or not wrapper._spoke):
                spoke = self.report(fmt % ret, level)
                if not wrapper._spoke:
                    wrapper._spoke = spoke
            return ret
        wrapper._spoke = False
        wrapper.__doc__ = func.__doc__
        return wrapper

    @cbook.deprecated("2.2", message=_verbose_msg)
    def ge(self, level):
        'return true if self.level is >= level'
        return self.vald[self.level] >= self.vald[level]


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    verbose = Verbose()


def _wrap(fmt, func, level=logging.DEBUG, always=True):
    """
    return a callable function that wraps func and reports its
    output through logger

    if always is True, the report will occur on every function
    call; otherwise only on the first time the function is called
    """
    assert callable(func)

    def wrapper(*args, **kwargs):
        ret = func(*args, **kwargs)

        if (always or not wrapper._spoke):
            _log.log(level, fmt % ret)
            spoke = True
            if not wrapper._spoke:
                wrapper._spoke = spoke
        return ret
    wrapper._spoke = False
    wrapper.__doc__ = func.__doc__
    return wrapper


def checkdep_dvipng():
    try:
        s = subprocess.Popen([str('dvipng'), '-version'],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        stdout, stderr = s.communicate()
        line = stdout.decode('ascii').split('\n')[1]
        v = line.split()[-1]
        return v
    except (IndexError, ValueError, OSError):
        return None


def checkdep_ghostscript():
    if checkdep_ghostscript.executable is None:
        if sys.platform == 'win32':
            # mgs is the name in miktex
            gs_execs = ['gswin32c', 'gswin64c', 'mgs', 'gs']
        else:
            gs_execs = ['gs']
        for gs_exec in gs_execs:
            try:
                s = subprocess.Popen(
                    [str(gs_exec), '--version'], stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
                stdout, stderr = s.communicate()
                if s.returncode == 0:
                    v = stdout[:-1].decode('ascii')
                    checkdep_ghostscript.executable = gs_exec
                    checkdep_ghostscript.version = v
            except (IndexError, ValueError, OSError):
                pass
    return checkdep_ghostscript.executable, checkdep_ghostscript.version
checkdep_ghostscript.executable = None
checkdep_ghostscript.version = None


# Deprecated, as it is unneeded and some distributions (e.g. MiKTeX 2.9.6350)
# do not actually report the TeX version.
@cbook.deprecated("2.1")
def checkdep_tex():
    try:
        s = subprocess.Popen([str('tex'), '-version'], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        stdout, stderr = s.communicate()
        line = stdout.decode('ascii').split('\n')[0]
        pattern = r'3\.1\d+'
        match = re.search(pattern, line)
        v = match.group(0)
        return v
    except (IndexError, ValueError, AttributeError, OSError):
        return None


def checkdep_pdftops():
    try:
        s = subprocess.Popen([str('pdftops'), '-v'], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        stdout, stderr = s.communicate()
        lines = stderr.decode('ascii').split('\n')
        for line in lines:
            if 'version' in line:
                v = line.split()[-1]
        return v
    except (IndexError, ValueError, UnboundLocalError, OSError):
        return None


def checkdep_inkscape():
    if checkdep_inkscape.version is None:
        try:
            s = subprocess.Popen([str('inkscape'), '-V'],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            stdout, stderr = s.communicate()
            lines = stdout.decode('ascii').split('\n')
            for line in lines:
                if 'Inkscape' in line:
                    v = line.split()[1]
                    break
            checkdep_inkscape.version = v
        except (IndexError, ValueError, UnboundLocalError, OSError):
            pass
    return checkdep_inkscape.version
checkdep_inkscape.version = None


@cbook.deprecated("2.1")
def checkdep_xmllint():
    try:
        s = subprocess.Popen([str('xmllint'), '--version'],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        stdout, stderr = s.communicate()
        lines = stderr.decode('ascii').split('\n')
        for line in lines:
            if 'version' in line:
                v = line.split()[-1]
                break
        return v
    except (IndexError, ValueError, UnboundLocalError, OSError):
        return None


def checkdep_ps_distiller(s):
    if not s:
        return False

    flag = True
    gs_req = '8.60'
    gs_exec, gs_v = checkdep_ghostscript()
    if not compare_versions(gs_v, gs_req):
        flag = False
        warnings.warn(('matplotlibrc ps.usedistiller option can not be used '
                       'unless ghostscript-%s or later is installed on your '
                       'system') % gs_req)

    if s == 'xpdf':
        pdftops_req = '3.0'
        pdftops_req_alt = '0.9'  # poppler version numbers, ugh
        pdftops_v = checkdep_pdftops()
        if compare_versions(pdftops_v, pdftops_req):
            pass
        elif (compare_versions(pdftops_v, pdftops_req_alt) and not
              compare_versions(pdftops_v, '1.0')):
            pass
        else:
            flag = False
            warnings.warn(('matplotlibrc ps.usedistiller can not be set to '
                           'xpdf unless xpdf-%s or later is installed on '
                           'your system') % pdftops_req)

    if flag:
        return s
    else:
        return False


def checkdep_usetex(s):
    if not s:
        return False

    gs_req = '8.60'
    dvipng_req = '1.6'
    flag = True

    if _backports.which("tex") is None:
        flag = False
        warnings.warn('matplotlibrc text.usetex option can not be used unless '
                      'TeX is installed on your system')

    dvipng_v = checkdep_dvipng()
    if not compare_versions(dvipng_v, dvipng_req):
        flag = False
        warnings.warn('matplotlibrc text.usetex can not be used with *Agg '
                      'backend unless dvipng-%s or later is installed on '
                      'your system' % dvipng_req)

    gs_exec, gs_v = checkdep_ghostscript()
    if not compare_versions(gs_v, gs_req):
        flag = False
        warnings.warn('matplotlibrc text.usetex can not be used unless '
                      'ghostscript-%s or later is installed on your system'
                      % gs_req)

    return flag


def _get_home():
    """Find user's home directory if possible.
    Otherwise, returns None.

    :see:
        http://mail.python.org/pipermail/python-list/2005-February/325395.html
    """
    if six.PY2 and sys.platform == 'win32':
        path = os.path.expanduser(b"~").decode(sys.getfilesystemencoding())
    else:
        path = os.path.expanduser("~")
    if os.path.isdir(path):
        return path
    for evar in ('HOME', 'USERPROFILE', 'TMP'):
        path = os.environ.get(evar)
        if path is not None and os.path.isdir(path):
            return path
    return None


def _create_tmp_config_dir():
    """
    If the config directory can not be created, create a temporary
    directory.
    """
    configdir = os.environ['MPLCONFIGDIR'] = (
        tempfile.mkdtemp(prefix='matplotlib-'))
    atexit.register(shutil.rmtree, configdir)
    return configdir


get_home = _wrap('$HOME=%s', _get_home, always=False)


def _get_xdg_config_dir():
    """
    Returns the XDG configuration directory, according to the `XDG
    base directory spec
    <http://standards.freedesktop.org/basedir-spec/basedir-spec-latest.html>`_.
    """
    path = os.environ.get('XDG_CONFIG_HOME')
    if path is None:
        path = get_home()
        if path is not None:
            path = os.path.join(path, '.config')
    return path


def _get_xdg_cache_dir():
    """
    Returns the XDG cache directory, according to the `XDG
    base directory spec
    <http://standards.freedesktop.org/basedir-spec/basedir-spec-latest.html>`_.
    """
    path = os.environ.get('XDG_CACHE_HOME')
    if path is None:
        path = get_home()
        if path is not None:
            path = os.path.join(path, '.cache')
    return path


def _get_config_or_cache_dir(xdg_base):
    from matplotlib.cbook import mkdirs

    configdir = os.environ.get('MPLCONFIGDIR')
    if configdir is not None:
        configdir = os.path.abspath(configdir)
        if not os.path.exists(configdir):
            mkdirs(configdir)

        if not _is_writable_dir(configdir):
            return _create_tmp_config_dir()
        return configdir

    p = None
    h = get_home()
    if h is not None:
        p = os.path.join(h, '.matplotlib')
    if sys.platform.startswith(('linux', 'freebsd')):
        p = None
        if xdg_base is not None:
            p = os.path.join(xdg_base, 'matplotlib')

    if p is not None:
        if os.path.exists(p):
            if _is_writable_dir(p):
                return p
        else:
            try:
                mkdirs(p)
            except OSError:
                pass
            else:
                return p

    return _create_tmp_config_dir()


def _get_configdir():
    """
    Return the string representing the configuration directory.

    The directory is chosen as follows:

    1. If the MPLCONFIGDIR environment variable is supplied, choose that.

    2a. On Linux, follow the XDG specification and look first in
        `$XDG_CONFIG_HOME`, if defined, or `$HOME/.config`.

    2b. On other platforms, choose `$HOME/.matplotlib`.

    3. If the chosen directory exists and is writable, use that as the
       configuration directory.
    4. If possible, create a temporary directory, and use it as the
       configuration directory.
    5. A writable directory could not be found or created; return None.
    """
    return _get_config_or_cache_dir(_get_xdg_config_dir())

get_configdir = _wrap('CONFIGDIR=%s', _get_configdir, always=False)


def _get_cachedir():
    """
    Return the location of the cache directory.

    The procedure used to find the directory is the same as for
    _get_config_dir, except using `$XDG_CACHE_HOME`/`~/.cache` instead.
    """
    return _get_config_or_cache_dir(_get_xdg_cache_dir())

get_cachedir = _wrap('CACHEDIR=%s', _get_cachedir, always=False)


def _decode_filesystem_path(path):
    if not isinstance(path, str):
        return path.decode(sys.getfilesystemencoding())
    else:
        return path


def _get_data_path():
    'get the path to matplotlib data'

    if 'MATPLOTLIBDATA' in os.environ:
        path = os.environ['MATPLOTLIBDATA']
        if not os.path.isdir(path):
            raise RuntimeError('Path in environment MATPLOTLIBDATA not a '
                               'directory')
        return path

    _file = _decode_filesystem_path(__file__)
    path = os.sep.join([os.path.dirname(_file), 'mpl-data'])
    if os.path.isdir(path):
        return path

    # setuptools' namespace_packages may highjack this init file
    # so need to try something known to be in matplotlib, not basemap
    import matplotlib.afm
    _file = _decode_filesystem_path(matplotlib.afm.__file__)
    path = os.sep.join([os.path.dirname(_file), 'mpl-data'])
    if os.path.isdir(path):
        return path

    # py2exe zips pure python, so still need special check
    if getattr(sys, 'frozen', None):
        exe_path = os.path.dirname(_decode_filesystem_path(sys.executable))
        path = os.path.join(exe_path, 'mpl-data')
        if os.path.isdir(path):
            return path

        # Try again assuming we need to step up one more directory
        path = os.path.join(os.path.split(exe_path)[0], 'mpl-data')
        if os.path.isdir(path):
            return path

        # Try again assuming sys.path[0] is a dir not a exe
        path = os.path.join(sys.path[0], 'mpl-data')
        if os.path.isdir(path):
            return path

    raise RuntimeError('Could not find the matplotlib data files')


def _get_data_path_cached():
    if defaultParams['datapath'][0] is None:
        defaultParams['datapath'][0] = _get_data_path()
    return defaultParams['datapath'][0]

get_data_path = _wrap('matplotlib data path %s', _get_data_path_cached,
                      always=False)


def get_py2exe_datafiles():
    datapath = get_data_path()
    _, tail = os.path.split(datapath)
    d = {}
    for root, _, files in os.walk(datapath):
        # Need to explicitly remove cocoa_agg files or py2exe complains
        # NOTE I don't know why, but do as previous version
        if 'Matplotlib.nib' in files:
            files.remove('Matplotlib.nib')
        files = [os.path.join(root, filename) for filename in files]
        root = root.replace(tail, 'mpl-data')
        root = root[root.index('mpl-data'):]
        d[root] = files
    return list(d.items())


def matplotlib_fname():
    """
    Get the location of the config file.

    The file location is determined in the following order

    - `$PWD/matplotlibrc`

    - `$MATPLOTLIBRC` if it is a file (or a named pipe, which can be created
      e.g. by process substitution)

    - `$MATPLOTLIBRC/matplotlibrc`

    - `$MPLCONFIGDIR/matplotlibrc`

    - On Linux,

          - `$XDG_CONFIG_HOME/matplotlib/matplotlibrc` (if
            $XDG_CONFIG_HOME is defined)

          - or `$HOME/.config/matplotlib/matplotlibrc` (if
            $XDG_CONFIG_HOME is not defined)

    - On other platforms,

         - `$HOME/.matplotlib/matplotlibrc` if `$HOME` is defined.

    - Lastly, it looks in `$MATPLOTLIBDATA/matplotlibrc` for a
      system-defined copy.
    """

    def gen_candidates():
        yield os.path.join(six.moves.getcwd(), 'matplotlibrc')
        try:
            matplotlibrc = os.environ['MATPLOTLIBRC']
        except KeyError:
            pass
        else:
            yield matplotlibrc
            yield os.path.join(matplotlibrc, 'matplotlibrc')
        yield os.path.join(_get_configdir(), 'matplotlibrc')
        yield os.path.join(get_data_path(), 'matplotlibrc')

    for fname in gen_candidates():
        if os.path.exists(fname):
            st_mode = os.stat(fname).st_mode
            if stat.S_ISREG(st_mode) or stat.S_ISFIFO(st_mode):
                break
    # Return first candidate that is a file, or last candidate if none is
    # valid (in that case, a warning is raised at startup by `rc_params`).
    return fname


# names of keys to deprecate
# the values are a tuple of (new_name, f_old_2_new, f_new_2_old)
# the inverse function may be `None`
_deprecated_map = {}

_deprecated_ignore_map = {'nbagg.transparent': 'figure.facecolor'}

_obsolete_set = {'plugins.directory', 'text.dvipnghack'}

# The following may use a value of None to suppress the warning.
# do NOT include in _all_deprecated
_deprecated_set = {'axes.hold',
                   'backend.qt4',
                   'backend.qt5'}

_all_deprecated = set(itertools.chain(
    _deprecated_ignore_map, _deprecated_map, _obsolete_set))


class RcParams(MutableMapping, dict):

    """
    A dictionary object including validation

    validating functions are defined and associated with rc parameters in
    :mod:`matplotlib.rcsetup`
    """

    validate = dict((key, converter) for key, (default, converter) in
                    six.iteritems(defaultParams)
                    if key not in _all_deprecated)
    msg_depr = "%s is deprecated and replaced with %s; please use the latter."
    msg_depr_set = ("%s is deprecated. Please remove it from your "
                    "matplotlibrc and/or style files.")
    msg_depr_ignore = "%s is deprecated and ignored. Use %s instead."
    msg_obsolete = ("%s is obsolete. Please remove it from your matplotlibrc "
                    "and/or style files.")
    msg_backend_obsolete = ("The {} rcParam was deprecated in version 2.2.  In"
                            " order to force the use of a specific Qt binding,"
                            " either import that binding first, or set the "
                            "QT_API environment variable.")

    # validate values on the way in
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __setitem__(self, key, val):
        try:
            if key in _deprecated_map:
                alt_key, alt_val, inverse_alt = _deprecated_map[key]
                warnings.warn(self.msg_depr % (key, alt_key),
                              mplDeprecation)
                key = alt_key
                val = alt_val(val)
            elif key in _deprecated_set and val is not None:
                if key.startswith('backend'):
                    warnings.warn(self.msg_backend_obsolete.format(key),
                                  mplDeprecation)
                else:
                    warnings.warn(self.msg_depr_set % key,
                                  mplDeprecation)
            elif key in _deprecated_ignore_map:
                alt = _deprecated_ignore_map[key]
                warnings.warn(self.msg_depr_ignore % (key, alt),
                              mplDeprecation)
                return
            elif key in _obsolete_set:
                warnings.warn(self.msg_obsolete % (key, ),
                              mplDeprecation)
                return
            try:
                cval = self.validate[key](val)
            except ValueError as ve:
                raise ValueError("Key %s: %s" % (key, str(ve)))
            dict.__setitem__(self, key, cval)
        except KeyError:
            raise KeyError(
                '%s is not a valid rc parameter. See rcParams.keys() for a '
                'list of valid parameters.' % (key,))

    def __getitem__(self, key):
        inverse_alt = None
        if key in _deprecated_map:
            alt_key, alt_val, inverse_alt = _deprecated_map[key]
            warnings.warn(self.msg_depr % (key, alt_key),
                          mplDeprecation)
            key = alt_key

        elif key in _deprecated_ignore_map:
            alt = _deprecated_ignore_map[key]
            warnings.warn(self.msg_depr_ignore % (key, alt),
                          mplDeprecation)
            key = alt

        elif key in _obsolete_set:
            warnings.warn(self.msg_obsolete % (key, ),
                          mplDeprecation)
            return None

        val = dict.__getitem__(self, key)
        if inverse_alt is not None:
            return inverse_alt(val)
        else:
            return val

    def __repr__(self):
        import pprint
        class_name = self.__class__.__name__
        indent = len(class_name) + 1
        repr_split = pprint.pformat(dict(self), indent=1,
                                    width=80 - indent).split('\n')
        repr_indented = ('\n' + ' ' * indent).join(repr_split)
        return '{0}({1})'.format(class_name, repr_indented)

    def __str__(self):
        return '\n'.join('{0}: {1}'.format(k, v)
                         for k, v in sorted(self.items()))

    def __iter__(self):
        """
        Yield sorted list of keys.
        """
        for k in sorted(dict.__iter__(self)):
            yield k

    def find_all(self, pattern):
        """
        Return the subset of this RcParams dictionary whose keys match,
        using :func:`re.search`, the given ``pattern``.

        .. note::

            Changes to the returned dictionary are *not* propagated to
            the parent RcParams dictionary.

        """
        pattern_re = re.compile(pattern)
        return RcParams((key, value)
                        for key, value in self.items()
                        if pattern_re.search(key))


def rc_params(fail_on_error=False):
    """Return a :class:`matplotlib.RcParams` instance from the
    default matplotlib rc file.
    """
    fname = matplotlib_fname()
    if not os.path.exists(fname):
        # this should never happen, default in mpl-data should always be found
        message = 'could not find rc file; returning defaults'
        ret = RcParams([(key, default) for key, (default, _) in
                        six.iteritems(defaultParams)
                        if key not in _all_deprecated])
        warnings.warn(message)
        return ret

    return rc_params_from_file(fname, fail_on_error)


URL_REGEX = re.compile(r'http://|https://|ftp://|file://|file:\\')


def is_url(filename):
    """Return True if string is an http, ftp, or file URL path."""
    return URL_REGEX.match(filename) is not None


def _url_lines(f):
    # Compatibility for urlopen in python 3, which yields bytes.
    for line in f:
        yield line.decode('utf8')


@contextlib.contextmanager
def _open_file_or_url(fname):
    if is_url(fname):
        f = urlopen(fname)
        yield _url_lines(f)
        f.close()
    else:
        fname = os.path.expanduser(fname)
        encoding = locale.getpreferredencoding(do_setlocale=False)
        if encoding is None:
            encoding = "utf-8"
        with io.open(fname, encoding=encoding) as f:
            yield f


_error_details_fmt = 'line #%d\n\t"%s"\n\tin file "%s"'


def _rc_params_in_file(fname, fail_on_error=False):
    """Return :class:`matplotlib.RcParams` from the contents of the given file.

    Unlike `rc_params_from_file`, the configuration class only contains the
    parameters specified in the file (i.e. default values are not filled in).
    """
    cnt = 0
    rc_temp = {}
    with _open_file_or_url(fname) as fd:
        try:
            for line in fd:
                cnt += 1
                strippedline = line.split('#', 1)[0].strip()
                if not strippedline:
                    continue
                tup = strippedline.split(':', 1)
                if len(tup) != 2:
                    error_details = _error_details_fmt % (cnt, line, fname)
                    warnings.warn('Illegal %s' % error_details)
                    continue
                key, val = tup
                key = key.strip()
                val = val.strip()
                if key in rc_temp:
                    warnings.warn('Duplicate key in file "%s", line #%d' %
                                  (fname, cnt))
                rc_temp[key] = (val, line, cnt)
        except UnicodeDecodeError:
            warnings.warn(
                ('Cannot decode configuration file %s with '
                 'encoding %s, check LANG and LC_* variables')
                % (fname, locale.getpreferredencoding(do_setlocale=False) or
                   'utf-8 (default)'))
            raise

    config = RcParams()

    for key in ('verbose.level', 'verbose.fileo'):
        if key in rc_temp:
            val, line, cnt = rc_temp.pop(key)
            if fail_on_error:
                config[key] = val  # try to convert to proper type or raise
            else:
                try:
                    config[key] = val  # try to convert to proper type or skip
                except Exception as msg:
                    error_details = _error_details_fmt % (cnt, line, fname)
                    warnings.warn('Bad val "%s" on %s\n\t%s' %
                                  (val, error_details, msg))

    for key, (val, line, cnt) in six.iteritems(rc_temp):
        if key in defaultParams:
            if fail_on_error:
                config[key] = val  # try to convert to proper type or raise
            else:
                try:
                    config[key] = val  # try to convert to proper type or skip
                except Exception as msg:
                    error_details = _error_details_fmt % (cnt, line, fname)
                    warnings.warn('Bad val "%s" on %s\n\t%s' %
                                  (val, error_details, msg))
        elif key in _deprecated_ignore_map:
            warnings.warn('%s is deprecated. Update your matplotlibrc to use '
                          '%s instead.' % (key, _deprecated_ignore_map[key]),
                          mplDeprecation)

        else:
            print("""
Bad key "%s" on line %d in
%s.
You probably need to get an updated matplotlibrc file from
http://github.com/matplotlib/matplotlib/blob/master/matplotlibrc.template
or from the matplotlib source distribution""" % (key, cnt, fname),
                  file=sys.stderr)

    return config


def rc_params_from_file(fname, fail_on_error=False, use_default_template=True):
    """Return :class:`matplotlib.RcParams` from the contents of the given file.

    Parameters
    ----------
    fname : str
        Name of file parsed for matplotlib settings.
    fail_on_error : bool
        If True, raise an error when the parser fails to convert a parameter.
    use_default_template : bool
        If True, initialize with default parameters before updating with those
        in the given file. If False, the configuration class only contains the
        parameters specified in the file. (Useful for updating dicts.)
    """
    config_from_file = _rc_params_in_file(fname, fail_on_error)

    if not use_default_template:
        return config_from_file

    iter_params = six.iteritems(defaultParams)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", mplDeprecation)
        config = RcParams([(key, default) for key, (default, _) in iter_params
                           if key not in _all_deprecated])
    config.update(config_from_file)

    if config['datapath'] is None:
        config['datapath'] = get_data_path()

    if "".join(config['text.latex.preamble']):
        _log.info("""
*****************************************************************
You have the following UNSUPPORTED LaTeX preamble customizations:
%s
Please do not ask for support with these customizations active.
*****************************************************************
""", '\n'.join(config['text.latex.preamble']))
    _log.debug('loaded rc file %s', fname)

    return config


# this is the instance used by the matplotlib classes
rcParams = rc_params()

if rcParams['examples.directory']:
    # paths that are intended to be relative to matplotlib_fname()
    # are allowed for the examples.directory parameter.
    # However, we will need to fully qualify the path because
    # Sphinx requires absolute paths.
    if not os.path.isabs(rcParams['examples.directory']):
        _basedir, _fname = os.path.split(matplotlib_fname())
        # Sometimes matplotlib_fname() can return relative paths,
        # Also, using realpath() guarantees that Sphinx will use
        # the same path that matplotlib sees (in case of weird symlinks).
        _basedir = os.path.realpath(_basedir)
        _fullpath = os.path.join(_basedir, rcParams['examples.directory'])
        rcParams['examples.directory'] = _fullpath

rcParamsOrig = rcParams.copy()

with warnings.catch_warnings():
    warnings.simplefilter("ignore", mplDeprecation)
    rcParamsDefault = RcParams([(key, default) for key, (default, converter) in
                                six.iteritems(defaultParams)
                                if key not in _all_deprecated])

rcParams['ps.usedistiller'] = checkdep_ps_distiller(
                      rcParams['ps.usedistiller'])

rcParams['text.usetex'] = checkdep_usetex(rcParams['text.usetex'])

if rcParams['axes.formatter.use_locale']:
    locale.setlocale(locale.LC_ALL, '')


def rc(group, **kwargs):
    """
    Set the current rc params.  Group is the grouping for the rc, e.g.,
    for ``lines.linewidth`` the group is ``lines``, for
    ``axes.facecolor``, the group is ``axes``, and so on.  Group may
    also be a list or tuple of group names, e.g., (*xtick*, *ytick*).
    *kwargs* is a dictionary attribute name/value pairs, e.g.,::

      rc('lines', linewidth=2, color='r')

    sets the current rc params and is equivalent to::

      rcParams['lines.linewidth'] = 2
      rcParams['lines.color'] = 'r'

    The following aliases are available to save typing for interactive
    users:

    =====   =================
    Alias   Property
    =====   =================
    'lw'    'linewidth'
    'ls'    'linestyle'
    'c'     'color'
    'fc'    'facecolor'
    'ec'    'edgecolor'
    'mew'   'markeredgewidth'
    'aa'    'antialiased'
    =====   =================

    Thus you could abbreviate the above rc command as::

          rc('lines', lw=2, c='r')


    Note you can use python's kwargs dictionary facility to store
    dictionaries of default parameters.  e.g., you can customize the
    font rc as follows::

      font = {'family' : 'monospace',
              'weight' : 'bold',
              'size'   : 'larger'}

      rc('font', **font)  # pass in the font dict as kwargs

    This enables you to easily switch between several configurations.  Use
    ``matplotlib.style.use('default')`` or :func:`~matplotlib.rcdefaults` to
    restore the default rc params after changes.
    """

    aliases = {
        'lw':  'linewidth',
        'ls':  'linestyle',
        'c':   'color',
        'fc':  'facecolor',
        'ec':  'edgecolor',
        'mew': 'markeredgewidth',
        'aa':  'antialiased',
        }

    if isinstance(group, six.string_types):
        group = (group,)
    for g in group:
        for k, v in six.iteritems(kwargs):
            name = aliases.get(k) or k
            key = '%s.%s' % (g, name)
            try:
                rcParams[key] = v
            except KeyError:
                raise KeyError(('Unrecognized key "%s" for group "%s" and '
                                'name "%s"') % (key, g, name))


def rcdefaults():
    """Restore the rc params from Matplotlib's internal defaults.

    See Also
    --------
    rc_file_defaults :
        Restore the rc params from the rc file originally loaded by Matplotlib.
    matplotlib.style.use :
        Use a specific style file.  Call ``style.use('default')`` to restore
        the default style.
    """
    rcParams.clear()
    rcParams.update(rcParamsDefault)


def rc_file_defaults():
    """Restore the rc params from the original rc file loaded by Matplotlib.
    """
    rcParams.update(rcParamsOrig)


def rc_file(fname):
    """
    Update rc params from file.
    """
    rcParams.update(rc_params_from_file(fname))


@contextlib.contextmanager
def rc_context(rc=None, fname=None):
    """
    Return a context manager for managing rc settings.

    This allows one to do::

        with mpl.rc_context(fname='screen.rc'):
            plt.plot(x, a)
            with mpl.rc_context(fname='print.rc'):
                plt.plot(x, b)
            plt.plot(x, c)

    The 'a' vs 'x' and 'c' vs 'x' plots would have settings from
    'screen.rc', while the 'b' vs 'x' plot would have settings from
    'print.rc'.

    A dictionary can also be passed to the context manager::

        with mpl.rc_context(rc={'text.usetex': True}, fname='screen.rc'):
            plt.plot(x, a)

    The 'rc' dictionary takes precedence over the settings loaded from
    'fname'.  Passing a dictionary only is also valid. For example a
    common usage is::

        with mpl.rc_context(rc={'interactive': False}):
            fig, ax = plt.subplots()
            ax.plot(range(3), range(3))
            fig.savefig('A.png', format='png')
            plt.close(fig)

    """

    orig = rcParams.copy()
    try:
        if fname:
            rc_file(fname)
        if rc:
            rcParams.update(rc)
        yield
    finally:
        # No need to revalidate the original values.
        dict.update(rcParams, orig)


_use_error_msg = """
This call to matplotlib.use() has no effect because the backend has already
been chosen; matplotlib.use() must be called *before* pylab, matplotlib.pyplot,
or matplotlib.backends is imported for the first time.

The backend was *originally* set to {backend!r} by the following code:
{tb}
"""


def use(arg, warn=True, force=False):
    """
    Set the matplotlib backend to one of the known backends.

    The argument is case-insensitive. *warn* specifies whether a
    warning should be issued if a backend has already been set up.
    *force* is an **experimental** flag that tells matplotlib to
    attempt to initialize a new backend by reloading the backend
    module.

    .. note::

        This function must be called *before* importing pyplot for
        the first time; or, if you are not using pyplot, it must be called
        before importing matplotlib.backends.  If warn is True, a warning
        is issued if you try and call this after pylab or pyplot have been
        loaded.  In certain black magic use cases, e.g.
        :func:`pyplot.switch_backend`, we are doing the reloading necessary to
        make the backend switch work (in some cases, e.g., pure image
        backends) so one can set warn=False to suppress the warnings.

    To find out which backend is currently set, see
    :func:`matplotlib.get_backend`.

    """
    # Lets determine the proper backend name first
    if arg.startswith('module://'):
        name = arg
    else:
        # Lowercase only non-module backend names (modules are case-sensitive)
        arg = arg.lower()
        name = validate_backend(arg)

    # Check if we've already set up a backend
    if 'matplotlib.backends' in sys.modules:
        # Warn only if called with a different name
        if (rcParams['backend'] != name) and warn:
            import matplotlib.backends
            warnings.warn(
                _use_error_msg.format(
                    backend=rcParams['backend'],
                    tb=matplotlib.backends._backend_loading_tb),
                stacklevel=2)

        # Unless we've been told to force it, just return
        if not force:
            return
        need_reload = True
    else:
        need_reload = False

    # Store the backend name
    rcParams['backend'] = name

    # If needed we reload here because a lot of setup code is triggered on
    # module import. See backends/__init__.py for more detail.
    if need_reload:
        reload(sys.modules['matplotlib.backends'])


try:
    use(os.environ['MPLBACKEND'])
except KeyError:
    pass


def get_backend():
    """Return the name of the current backend."""
    return rcParams['backend']


def interactive(b):
    """
    Set interactive mode to boolean b.

    If b is True, then draw after every plotting command, e.g., after xlabel
    """
    rcParams['interactive'] = b


def is_interactive():
    'Return true if plot mode is interactive'
    return rcParams['interactive']


def tk_window_focus():
    """Return true if focus maintenance under TkAgg on win32 is on.
     This currently works only for python.exe and IPython.exe.
     Both IDLE and Pythonwin.exe fail badly when tk_window_focus is on."""
    if rcParams['backend'] != 'TkAgg':
        return False
    return rcParams['tk.window_focus']


default_test_modules = [
    'matplotlib.tests',
    'matplotlib.sphinxext.tests',
    'mpl_toolkits.tests',
]


def _init_tests():
    try:
        import faulthandler
    except ImportError:
        pass
    else:
        # CPython's faulthandler since v3.6 handles exceptions on Windows
        # https://bugs.python.org/issue23848 but until v3.6.4 it was
        # printing non-fatal exceptions https://bugs.python.org/issue30557
        import platform
        if not (sys.platform == 'win32' and
                (3, 6) < sys.version_info < (3, 6, 4) and
                platform.python_implementation() == 'CPython'):
            faulthandler.enable()

    # The version of FreeType to install locally for running the
    # tests.  This must match the value in `setupext.py`
    LOCAL_FREETYPE_VERSION = '2.6.1'

    from matplotlib import ft2font
    if (ft2font.__freetype_version__ != LOCAL_FREETYPE_VERSION or
        ft2font.__freetype_build_type__ != 'local'):
        warnings.warn(
            "Matplotlib is not built with the correct FreeType version to run "
            "tests.  Set local_freetype=True in setup.cfg and rebuild. "
            "Expect many image comparison failures below. "
            "Expected freetype version {0}. "
            "Found freetype version {1}. "
            "Freetype build type is {2}local".format(
                LOCAL_FREETYPE_VERSION,
                ft2font.__freetype_version__,
                "" if ft2font.__freetype_build_type__ == 'local' else "not "
            )
        )

    try:
        import pytest
        try:
            from unittest import mock
        except ImportError:
            import mock
    except ImportError:
        print("matplotlib.test requires pytest and mock to run.")
        raise


def test(verbosity=None, coverage=False, switch_backend_warn=True,
         recursionlimit=0, **kwargs):
    """run the matplotlib test suite"""
    _init_tests()
    if not os.path.isdir(os.path.join(os.path.dirname(__file__), 'tests')):
        raise ImportError("Matplotlib test data is not installed")

    old_backend = get_backend()
    old_recursionlimit = sys.getrecursionlimit()
    try:
        use('agg')
        if recursionlimit:
            sys.setrecursionlimit(recursionlimit)
        import pytest

        args = kwargs.pop('argv', [])
        provide_default_modules = True
        use_pyargs = True
        for arg in args:
            if any(arg.startswith(module_path)
                   for module_path in default_test_modules):
                provide_default_modules = False
                break
            if os.path.exists(arg):
                provide_default_modules = False
                use_pyargs = False
                break
        if use_pyargs:
            args += ['--pyargs']
        if provide_default_modules:
            args += default_test_modules

        if coverage:
            args += ['--cov']

        if verbosity:
            args += ['-' + 'v' * verbosity]

        retcode = pytest.main(args, **kwargs)
    finally:
        if old_backend.lower() != 'agg':
            use(old_backend, warn=switch_backend_warn)
        if recursionlimit:
            sys.setrecursionlimit(old_recursionlimit)

    return retcode


test.__test__ = False  # pytest: this function is not a test


def _replacer(data, key):
    """Either returns data[key] or passes data back. Also
    converts input data to a sequence as needed.
    """
    # if key isn't a string don't bother
    if not isinstance(key, six.string_types):
        return (key)
    # try to use __getitem__
    try:
        return sanitize_sequence(data[key])
    # key does not exist, silently fall back to key
    except KeyError:
        return key


_DATA_DOC_APPENDIX = """

.. note::
    In addition to the above described arguments, this function can take a
    **data** keyword argument. If such a **data** argument is given, the
    following arguments are replaced by **data[<arg>]**:

    {replaced}
"""


def _add_data_doc(docstring, replace_names, replace_all_args):
    """Add documentation for a *data* field to the given docstring.

    Parameters
    ----------
    docstring : str
        The input docstring.
    replace_names : list of strings or None
        The list of parameter names which arguments should be replaced by
        `data[name]`. If None, all arguments are replaced if they are
        included in `data`.
    replace_all_args : bool
        If True, all arguments in *args get replaced, even if they are not
        in replace_names.

    Returns
    -------
        The augmented docstring.
    """
    if docstring is None:
        docstring = ''
    else:
        docstring = dedent(docstring)
    _repl = ""
    if replace_names is None:
        _repl = "* All positional and all keyword arguments."
    else:
        if len(replace_names) != 0:
            _repl = "* All arguments with the following names: '{names}'."
        if replace_all_args:
            _repl += "\n    * All positional arguments."
        _repl = _repl.format(names="', '".join(sorted(replace_names)))
    return docstring + _DATA_DOC_APPENDIX.format(replaced=_repl)


def _preprocess_data(replace_names=None, replace_all_args=False,
                     label_namer=None, positional_parameter_names=None):
    """
    A decorator to add a 'data' kwarg to any a function.  The signature
    of the input function must include the ax argument at the first position ::

       def foo(ax, *args, **kwargs)

    so this is suitable for use with Axes methods.

    Parameters
    ----------
    replace_names : list of strings, optional, default: None
        The list of parameter names which arguments should be replaced by
        `data[name]`. If None, all arguments are replaced if they are
        included in `data`.
    replace_all_args : bool, default: False
        If True, all arguments in *args get replaced, even if they are not
        in replace_names.
    label_namer : string, optional, default: None
        The name of the parameter which argument should be used as label, if
        label is not set. If None, the label keyword argument is not set.
    positional_parameter_names : list of strings or callable, optional
        The full list of positional parameter names (excluding an explicit
        `ax`/'self' argument at the first place and including all possible
        positional parameter in `*args`), in the right order. Can also include
        all other keyword parameter. Only needed if the wrapped function does
        contain `*args` and (replace_names is not None or replace_all_args is
        False). If it is a callable, it will be called with the actual
        tuple of *args and the data and should return a list like
        above.
        NOTE: callables should only be used when the names and order of *args
        can only be determined at runtime. Please use list of names
        when the order and names of *args is clear before runtime!

    .. note:: decorator also converts MappingView input data to list.
    """
    if replace_names is not None:
        replace_names = set(replace_names)

    def param(func):
        new_sig = None
        # signature is since 3.3 and wrapped since 3.2, but we support 3.4+.
        python_has_signature = python_has_wrapped = six.PY3

        # if in a legacy version of python and IPython is already imported
        # try to use their back-ported signature
        if not python_has_signature and 'IPython' in sys.modules:
            try:
                import IPython.utils.signatures
                signature = IPython.utils.signatures.signature
                Parameter = IPython.utils.signatures.Parameter
            except ImportError:
                pass
            else:
                python_has_signature = True
        else:
            if python_has_signature:
                signature = inspect.signature
                Parameter = inspect.Parameter

        if not python_has_signature:
            arg_spec = inspect.getargspec(func)
            _arg_names = arg_spec.args
            _has_varargs = arg_spec.varargs is not None
            _has_varkwargs = arg_spec.keywords is not None
        else:
            sig = signature(func)
            _has_varargs = False
            _has_varkwargs = False
            _arg_names = []
            params = list(sig.parameters.values())
            for p in params:
                if p.kind is Parameter.VAR_POSITIONAL:
                    _has_varargs = True
                elif p.kind is Parameter.VAR_KEYWORD:
                    _has_varkwargs = True
                else:
                    _arg_names.append(p.name)
            data_param = Parameter('data',
                                   Parameter.KEYWORD_ONLY,
                                   default=None)
            if _has_varkwargs:
                params.insert(-1, data_param)
            else:
                params.append(data_param)
            new_sig = sig.replace(parameters=params)
        # Import-time check: do we have enough information to replace *args?
        arg_names_at_runtime = False
        # there can't be any positional arguments behind *args and no
        # positional args can end up in **kwargs, so only *varargs make
        # problems.
        # http://stupidpythonideas.blogspot.de/2013/08/arguments-and-parameters.html
        if not _has_varargs:
            # all args are "named", so no problem
            # remove the first "ax" / self arg
            arg_names = _arg_names[1:]
        else:
            # Here we have "unnamed" variables and we need a way to determine
            # whether to replace a arg or not
            if replace_names is None:
                # all argnames should be replaced
                arg_names = None
            elif len(replace_names) == 0:
                # No argnames should be replaced
                arg_names = []
            elif len(_arg_names) > 1 and (positional_parameter_names is None):
                # we got no manual parameter names but more than an 'ax' ...
                if len(replace_names - set(_arg_names[1:])) == 0:
                    # all to be replaced arguments are in the list
                    arg_names = _arg_names[1:]
                else:
                    raise AssertionError(
                        "Got unknown 'replace_names' and wrapped function "
                        "{!r} uses '*args', need 'positional_parameter_names'"
                        .format(func.__name__))
            else:
                if positional_parameter_names is not None:
                    if callable(positional_parameter_names):
                        # determined by the function at runtime
                        arg_names_at_runtime = True
                        # so that we don't compute the label_pos at import time
                        arg_names = []
                    else:
                        arg_names = positional_parameter_names
                else:
                    if replace_all_args:
                        arg_names = []
                    else:
                        raise AssertionError(
                            "Got 'replace_names' and wrapped function {!r} "
                            "uses *args, need 'positional_parameter_names' or "
                            "'replace_all_args'".format(func.__name__))

        # compute the possible label_namer and label position in positional
        # arguments
        label_pos = 9999  # bigger than all "possible" argument lists
        label_namer_pos = 9999  # bigger than all "possible" argument lists
        if (label_namer and  # we actually want a label here ...
                arg_names and  # and we can determine a label in *args ...
                (label_namer in arg_names)):  # and it is in *args
            label_namer_pos = arg_names.index(label_namer)
            if "label" in arg_names:
                label_pos = arg_names.index("label")

        # Check the case we know a label_namer but we can't find it the
        # arg_names... Unfortunately the label_namer can be in **kwargs,
        # which we can't detect here and which results in a non-set label
        # which might surprise the user :-(
        if label_namer and not arg_names_at_runtime and not _has_varkwargs:
            if not arg_names:
                raise AssertionError(
                    "label_namer {!r} can't be found as the parameter without "
                    "'positional_parameter_names'".format(label_namer))
            elif label_namer not in arg_names:
                raise AssertionError(
                    "label_namer {!r} can't be found in the parameter names "
                    "(known argnames: %s).".format(label_namer, arg_names))
            else:
                # this is the case when the name is in arg_names
                pass

        @functools.wraps(func)
        def inner(ax, *args, **kwargs):
            # this is needed because we want to change these values if
            # arg_names_at_runtime==True, but python does not allow assigning
            # to a variable in a outer scope. So use some new local ones and
            # set them to the already computed values.
            _label_pos = label_pos
            _label_namer_pos = label_namer_pos
            _arg_names = arg_names

            label = None

            data = kwargs.pop('data', None)

            if data is None:  # data validation
                args = tuple(sanitize_sequence(a) for a in args)
            else:
                if arg_names_at_runtime:
                    # update the information about replace names and
                    # label position
                    _arg_names = positional_parameter_names(args, data)
                    if (label_namer and  # we actually want a label here ...
                            _arg_names and  # and we can find a label in *args
                            (label_namer in _arg_names)):  # and it is in *args
                        _label_namer_pos = _arg_names.index(label_namer)
                        if "label" in _arg_names:
                            _label_pos = arg_names.index("label")

                # save the current label_namer value so that it can be used as
                # a label
                if _label_namer_pos < len(args):
                    label = args[_label_namer_pos]
                else:
                    label = kwargs.get(label_namer, None)
                # ensure a string, as label can't be anything else
                if not isinstance(label, six.string_types):
                    label = None

                if (replace_names is None) or (replace_all_args is True):
                    # all should be replaced
                    args = tuple(_replacer(data, a) for
                                 j, a in enumerate(args))
                else:
                    # An arg is replaced if the arg_name of that position is
                    #   in replace_names ...
                    if len(_arg_names) < len(args):
                        raise RuntimeError(
                            "Got more args than function expects")
                    args = tuple(_replacer(data, a)
                                 if _arg_names[j] in replace_names else a
                                 for j, a in enumerate(args))

                if replace_names is None:
                    # replace all kwargs ...
                    kwargs = dict((k, _replacer(data, v))
                                  for k, v in six.iteritems(kwargs))
                else:
                    # ... or only if a kwarg of that name is in replace_names
                    kwargs = dict((k, _replacer(data, v)
                                   if k in replace_names else v)
                                  for k, v in six.iteritems(kwargs))

            # replace the label if this func "wants" a label arg and the user
            # didn't set one. Note: if the user puts in "label=None", it does
            # *NOT* get replaced!
            user_supplied_label = (
                (len(args) >= _label_pos) or  # label is included in args
                ('label' in kwargs)  # ... or in kwargs
            )
            if (label_namer and not user_supplied_label):
                if _label_namer_pos < len(args):
                    kwargs['label'] = get_label(args[_label_namer_pos], label)
                elif label_namer in kwargs:
                    kwargs['label'] = get_label(kwargs[label_namer], label)
                else:
                    warnings.warn(
                        "Tried to set a label via parameter %r in func %r but "
                        "couldn't find such an argument.\n"
                        "(This is a programming error, please report to "
                        "the Matplotlib list!)" % (label_namer, func.__name__),
                        RuntimeWarning, stacklevel=2)
            return func(ax, *args, **kwargs)

        inner.__doc__ = _add_data_doc(inner.__doc__,
                                      replace_names, replace_all_args)
        if not python_has_wrapped:
            inner.__wrapped__ = func
        if new_sig is not None:
            inner.__signature__ = new_sig
        return inner

    return param

_log.debug('matplotlib version %s', __version__)
_log.debug('interactive is %s', is_interactive())
_log.debug('platform is %s', sys.platform)
_log.debug('loaded modules: %s', list(sys.modules))
