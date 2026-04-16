"""OpenMP wrapper using a libgomp dynamically loaded library."""

from ctypes.util import find_library
from subprocess import check_output, CalledProcessError, DEVNULL
import ctypes
import os
import sys
import sysconfig

try:
    # there may be an environ modification when loading config
    from pythran.config import compiler
except ImportError:
    def compiler():
        return os.environ.get('CXX', 'c++')
cxx = compiler()


# This function and the `msvc_runtime_*` ones below are taken over from
# numpy.distutils
def get_shared_lib_extension(is_python_ext=False):
    """Return the correct file extension for shared libraries.

    Parameters
    ----------
    is_python_ext : bool, optional
        Whether the shared library is a Python extension.  Default is False.

    Returns
    -------
    so_ext : str
        The shared library extension.

    Notes
    -----
    For Python shared libs, `so_ext` will typically be '.so' on Linux and OS X,
    and '.pyd' on Windows.  For Python >= 3.2 `so_ext` has a tag prepended on
    POSIX systems according to PEP 3149.

    """
    confvars = sysconfig.get_config_vars()
    so_ext = confvars.get('EXT_SUFFIX', '')

    if not is_python_ext:
        # hardcode some known values to avoid some old distutils bugs
        if (sys.platform.startswith('linux') or
            sys.platform.startswith('gnukfreebsd')):
            so_ext = '.so'
        elif sys.platform.startswith('darwin'):
            so_ext = '.dylib'
        elif sys.platform.startswith('win'):
            so_ext = '.dll'
        else:
            # don't use long extension (e.g., .cpython-310-x64_64-linux-gnu.so',
            # see PEP 3149), but subtract that Python-specific part here
            so_ext = so_ext.replace('.' + confvars.get('SOABI'), '', 1)

    return so_ext


def msvc_runtime_version():
    "Return version of MSVC runtime library, as defined by __MSC_VER__ macro"
    msc_pos = sys.version.find('MSC v.')
    if msc_pos != -1:
        msc_ver = int(sys.version[msc_pos+6:msc_pos+10])
    else:
        msc_ver = None
    return msc_ver


def msvc_runtime_major():
    "Return major version of MSVC runtime coded like get_build_msvc_version"
    major = {1300:  70,  # MSVC 7.0
             1310:  71,  # MSVC 7.1
             1400:  80,  # MSVC 8
             1500:  90,  # MSVC 9  (aka 2008)
             1600: 100,  # MSVC 10 (aka 2010)
             1900: 140,  # MSVC 14 (aka 2015)
             1910: 141,  # MSVC 141 (aka 2017)
             1920: 142,  # MSVC 142 (aka 2019)
    }.get(msvc_runtime_version(), None)
    return major


class OpenMP(object):

    """
    Internal representation of the OpenMP module.

    Custom class is used to dynamically add omp runtime function
    to this library when function is called.
    """

    def __init__(self):
        # FIXME: this is broken, for at least two reasons:
        #   (1) checking how Python was built is not a good way to determine if
        #       MSVC is being used right now,
        #   (2) the `msvc_runtime_major` function above came from
        #       `numpy.distutils`, where it did not have entries for 1910/1920
        #       and returned None instead
        ver = msvc_runtime_major()
        if ver is None:
            self.init_not_msvc()
        else:
            self.init_msvc(ver)

    def init_msvc(self, ver):
        vcomp_path = find_library('vcomp%d.dll' % ver)
        if not vcomp_path:
            raise ImportError("I can't find a shared library for vcomp.")
        else:
            # Load the library (shouldn't fail with an absolute path right?)
            self.libomp = ctypes.CDLL(vcomp_path)
            self.version = 20

    def get_libomp_names(self):
        """Return list of OpenMP libraries to try"""
        return ['omp', 'gomp', 'iomp5']

    def init_not_msvc(self):
        """ Find OpenMP library and try to load if using ctype interface. """
        # find_library() does not automatically search LD_LIBRARY_PATH
        # until Python 3.6+, so we explicitly add it.
        # LD_LIBRARY_PATH is used on Linux, while macOS uses DYLD_LIBRARY_PATH
        # and DYLD_FALLBACK_LIBRARY_PATH.
        env_vars = []
        if sys.platform == 'darwin':
            env_vars = ['DYLD_LIBRARY_PATH', 'DYLD_FALLBACK_LIBRARY_PATH']
        else:
            env_vars = ['LD_LIBRARY_PATH']

        paths = []
        for env_var in env_vars:
            env_paths = os.environ.get(env_var, '')
            if env_paths:
                paths.extend(env_paths.split(os.pathsep))


        libomp_names = self.get_libomp_names()

        if cxx is not None:
            for libomp_name in libomp_names:
                cmd = [cxx,
                       '-print-file-name=lib{}{}'.format(
                           libomp_name,
                           get_shared_lib_extension())]
                # The subprocess can fail in various ways, including because it
                # doesn't support '-print-file-name'. In that case just give up.
                try:
                    output = check_output(cmd,
                                          stderr=DEVNULL)
                    path = os.path.dirname(output.decode().strip())
                    if path:
                        paths.append(path)
                except (OSError, CalledProcessError):
                    pass


        for libomp_name in libomp_names:
            # Try to load find libomp shared library using loader search dirs
            libomp_path = find_library(libomp_name)

            # Try to use custom paths if lookup failed
            if not libomp_path:
                for path in paths:
                    candidate_path = os.path.join(
                        path,
                        'lib{}{}'.format(libomp_name,
                                         get_shared_lib_extension()))
                    if os.path.isfile(candidate_path):
                        libomp_path = candidate_path
                        break

            # Load the library
            if libomp_path:
                try:
                    self.libomp = ctypes.CDLL(libomp_path)
                except OSError:
                    raise ImportError("found openMP library '{}' but couldn't load it. "
                                      "This may happen if you are cross-compiling.".format(libomp_path))
                self.version = 45
                return

        raise ImportError("I can't find a shared library for libomp, you may need to install it "
                          "or adjust the {} environment variable.".format(env_vars[0]))


    def __getattr__(self, name):
        """
        Get correct function name from libgomp ready to be use.

        __getattr__ is call only `name != libomp` as libomp is a real
        attribute.
        """
        if name == 'VERSION':
            return self.version
        return getattr(self.libomp, 'omp_' + name)

# see http://mail.python.org/pipermail/python-ideas/2012-May/014969.html
sys.modules[__name__] = OpenMP()

