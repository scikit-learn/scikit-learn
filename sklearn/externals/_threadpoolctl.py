"""threadpoolctl

This module provides utilities to introspect native libraries that relies on
thread pools (notably BLAS and OpenMP implementations) and dynamically set the
maximal number of threads they can use.
"""
# License: BSD 3-Clause

# The code to introspect dynamically loaded libraries on POSIX systems is
# adapted from code by Intel developper @anton-malakhov available at
# https://github.com/IntelPython/smp (Copyright (c) 2017, Intel Corporation)
# and also published under the BSD 3-Clause license
import os
import re
import sys
import ctypes
import textwrap
import warnings
from ctypes.util import find_library
from abc import ABC, abstractmethod

__version__ = "2.0.0"
__all__ = ["threadpool_limits", "threadpool_info"]


# One can get runtime errors or even segfaults due to multiple OpenMP libraries
# loaded simultaneously which can happen easily in Python when importing and
# using compiled extensions built with different compilers and therefore
# different OpenMP runtimes in the same program. In particular libiomp (used by
# Intel ICC) and libomp used by clang/llvm tend to crash. This can happen for
# instance when calling BLAS inside a prange. Setting the following environment
# variable allows multiple OpenMP libraries to be loaded. It should not degrade
# performances since we manually take care of potential over-subscription
# performance issues, in sections of the code where nested OpenMP loops can
# happen, by dynamically reconfiguring the inner OpenMP runtime to temporarily
# disable it while under the scope of the outer OpenMP parallel section.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "True")

# Structure to cast the info on dynamically loaded library. See
# https://linux.die.net/man/3/dl_iterate_phdr for more details.
_SYSTEM_UINT = ctypes.c_uint64 if sys.maxsize > 2**32 else ctypes.c_uint32
_SYSTEM_UINT_HALF = ctypes.c_uint32 if sys.maxsize > 2**32 else ctypes.c_uint16


class _dl_phdr_info(ctypes.Structure):
    _fields_ = [
        ("dlpi_addr",  _SYSTEM_UINT),       # Base address of object
        ("dlpi_name",  ctypes.c_char_p),    # path to the library
        ("dlpi_phdr",  ctypes.c_void_p),    # pointer on dlpi_headers
        ("dlpi_phnum",  _SYSTEM_UINT_HALF)  # number of elements in dlpi_phdr
    ]


# The RTLD_NOLOAD flag for loading shared libraries is not defined on Windows.
try:
    _RTLD_NOLOAD = os.RTLD_NOLOAD
except AttributeError:
    _RTLD_NOLOAD = ctypes.DEFAULT_MODE


# List of the supported libraries. The items are indexed by the name of the
# class to instanciate to create the module objects. The items hold the
# possible prefixes of loaded shared objects, the name of the internal_api to
# call and the name of the user_api.
_SUPPORTED_MODULES = {
    "_OpenMPModule": {
        "user_api": "openmp",
        "internal_api": "openmp",
        "filename_prefixes": ("libiomp", "libgomp", "libomp", "vcomp")
    },
    "_OpenBLASModule": {
        "user_api": "blas",
        "internal_api": "openblas",
        "filename_prefixes": ("libopenblas",)
    },
    "_MKLModule": {
        "user_api": "blas",
        "internal_api": "mkl",
        "filename_prefixes": ("libmkl_rt", "mkl_rt")
    },
    "_BLISModule": {
        "user_api": "blas",
        "internal_api": "blis",
        "filename_prefixes": ("libblis",)
    }
}

# Helpers for the doc and test names
_ALL_USER_APIS = list(set(m["user_api"] for m in _SUPPORTED_MODULES.values()))
_ALL_INTERNAL_APIS = [m["internal_api"] for m in _SUPPORTED_MODULES.values()]
_ALL_PREFIXES = [prefix for m in _SUPPORTED_MODULES.values()
                 for prefix in m["filename_prefixes"]]
_ALL_BLAS_LIBRARIES = [m["internal_api"] for m in _SUPPORTED_MODULES.values()
                       if m["user_api"] == "blas"]
_ALL_OPENMP_LIBRARIES = list(
    _SUPPORTED_MODULES["_OpenMPModule"]["filename_prefixes"])


def _format_docstring(*args, **kwargs):
    def decorator(o):
        o.__doc__ = o.__doc__.format(*args, **kwargs)
        return o

    return decorator


@_format_docstring(USER_APIS=list(_ALL_USER_APIS),
                   INTERNAL_APIS=_ALL_INTERNAL_APIS)
def threadpool_info():
    """Return the maximal number of threads for each detected library.

    Return a list with all the supported modules that have been found. Each
    module is represented by a dict with the following information:

      - "user_api" : user API. Possible values are {USER_APIS}.
      - "internal_api": internal API. Possible values are {INTERNAL_APIS}.
      - "prefix" : filename prefix of the specific implementation.
      - "filepath": path to the loaded module.
      - "version": version of the library (if available).
      - "num_threads": the current thread limit.

    In addition, each module may contain internal_api specific entries.
    """
    return _ThreadpoolInfo(user_api=_ALL_USER_APIS).todicts()


@_format_docstring(
    USER_APIS=", ".join('"{}"'.format(api) for api in _ALL_USER_APIS),
    BLAS_LIBS=", ".join(_ALL_BLAS_LIBRARIES),
    OPENMP_LIBS=", ".join(_ALL_OPENMP_LIBRARIES))
class threadpool_limits:
    """Change the maximal number of threads that can be used in thread pools.

    This class can be used either as a function (the construction of this
    object limits the number of threads) or as a context manager, in a `with`
    block.

    Set the maximal number of threads that can be used in thread pools used in
    the supported libraries to `limit`. This function works for libraries that
    are already loaded in the interpreter and can be changed dynamically.

    Parameters
    ----------
    limits : int, dict or None (default=None)
        The maximal number of threads that can be used in thread pools

        - If int, sets the maximum number of threads to `limits` for each
          library selected by `user_api`.

        - If it is a dictionary `{{key: max_threads}}`, this function sets a
          custom maximum number of threads for each `key` which can be either a
          `user_api` or a `prefix` for a specific library.

        - If None, this function does not do anything.

    user_api : {USER_APIS} or None (default=None)
        APIs of libraries to limit. Used only if `limits` is an int.

        - If "blas", it will only limit BLAS supported libraries ({BLAS_LIBS}).

        - If "openmp", it will only limit OpenMP supported libraries
          ({OPENMP_LIBS}). Note that it can affect the number of threads used
          by the BLAS libraries if they rely on OpenMP.

        - If None, this function will apply to all supported libraries.
    """
    def __init__(self, limits=None, user_api=None):
        self._limits, self._user_api, self._prefixes = \
            self._check_params(limits, user_api)

        self._original_info = self._set_threadpool_limits()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.unregister()

    def unregister(self):
        if self._original_info is not None:
            for module in self._original_info:
                module.set_num_threads(module.num_threads)

    def get_original_num_threads(self):
        """Original num_threads from before calling threadpool_limits

        Return a dict `{user_api: num_threads}`.
        """
        if self._original_info is not None:
            original_limits = self._original_info
        else:
            original_limits = _ThreadpoolInfo(user_api=self._user_api)

        num_threads = {}
        warning_apis = []

        for user_api in self._user_api:
            limits = [module.num_threads for module in
                      original_limits.get_modules("user_api", user_api)]
            limits = set(limits)
            n_limits = len(limits)

            if n_limits == 1:
                limit = limits.pop()
            elif n_limits == 0:
                limit = None
            else:
                limit = min(limits)
                warning_apis.append(user_api)

            num_threads[user_api] = limit

        if warning_apis:
            warnings.warn(
                "Multiple value possible for following user apis: "
                + ", ".join(warning_apis) + ". Returning the minimum.")

        return num_threads

    def _check_params(self, limits, user_api):
        """Suitable values for the _limits, _user_api and _prefixes attributes
        """
        if limits is None or isinstance(limits, int):
            if user_api is None:
                user_api = _ALL_USER_APIS
            elif user_api in _ALL_USER_APIS:
                user_api = [user_api]
            else:
                raise ValueError(
                    "user_api must be either in {} or None. Got "
                    "{} instead.".format(_ALL_USER_APIS, user_api))

            if limits is not None:
                limits = {api: limits for api in user_api}
            prefixes = []
        else:
            if isinstance(limits, list):
                # This should be a list of dicts of modules, for compatibility
                # with the result from threadpool_info.
                limits = {module["prefix"]: module["num_threads"]
                          for module in limits}
            elif isinstance(limits, _ThreadpoolInfo):
                # To set the limits from the modules of a _ThreadpoolInfo
                # object.
                limits = {module.prefix: module.num_threads
                          for module in limits}

            if not isinstance(limits, dict):
                raise TypeError("limits must either be an int, a list or a "
                                "dict. Got {} instead".format(type(limits)))

            # With a dictionary, can set both specific limit for given modules
            # and global limit for user_api. Fetch each separately.
            prefixes = [prefix for prefix in limits if prefix in _ALL_PREFIXES]
            user_api = [api for api in limits if api in _ALL_USER_APIS]

        return limits, user_api, prefixes

    def _set_threadpool_limits(self):
        """Change the maximal number of threads in selected thread pools.

        Return a list with all the supported modules that have been found
        matching `self._prefixes` and `self._user_api`.
        """
        if self._limits is None:
            return None

        modules = _ThreadpoolInfo(prefixes=self._prefixes,
                                  user_api=self._user_api)
        for module in modules:
            # self._limits is a dict {key: num_threads} where key is either
            # a prefix or a user_api. If a module matches both, the limit
            # corresponding to the prefix is chosed.
            if module.prefix in self._limits:
                num_threads = self._limits[module.prefix]
            else:
                num_threads = self._limits[module.user_api]

            if num_threads is not None:
                module.set_num_threads(num_threads)
        return modules


# The object oriented API of _ThreadpoolInfo and its modules is private.
# The public API (i.e. the "threadpool_info" function) only exposes the
# "list of dicts" representation returned by the .todicts method.
@_format_docstring(
    PREFIXES=", ".join('"{}"'.format(prefix) for prefix in _ALL_PREFIXES),
    USER_APIS=", ".join('"{}"'.format(api) for api in _ALL_USER_APIS),
    BLAS_LIBS=", ".join(_ALL_BLAS_LIBRARIES),
    OPENMP_LIBS=", ".join(_ALL_OPENMP_LIBRARIES))
class _ThreadpoolInfo():
    """Collection of all supported modules that have been found

    Parameters
    ----------
    user_api : list of user APIs or None (default=None)
        Select libraries matching the requested API. Ignored if `modules` is
        not None. Supported user APIs are {USER_APIS}.

        - "blas" selects all BLAS supported libraries ({BLAS_LIBS})
        - "openmp" selects all OpenMP supported libraries ({OPENMP_LIBS})

        If None, libraries are not selected by their `user_api`.

    prefixes : list of prefixes or None (default=None)
        Select libraries matching the requested prefixes. Supported prefixes
        are {PREFIXES}.
        If None, libraries are not selected by their prefix. Ignored if
        `modules` is not None.

    modules : list of _Module objects or None (default=None)
        Wraps a list of _Module objects into a _ThreapoolInfo object. Does not
        load or reload any shared library. If it is not None, `prefixes` and
        `user_api` are ignored.

    Note
    ----
    Is is possible to select libraries both by prefixes and by user_api. All
    libraries matching one or the other will be selected.
    """
    # Cache for libc under POSIX and a few system libraries under Windows.
    # We use a class level cache instead of an instance level cache because
    # it's very unlikely that a shared library will be unloaded and reloaded
    # during the lifetime of a program.
    _system_libraries = dict()
    # Cache for calls to os.path.realpath on system libraries to reduce the
    # impact of slow system calls (e.g. stat) on slow filesystem.
    # We use a class level cache instead of an instance level cache because
    # we can safely assume that the filepath of loaded shared libraries will
    # never change during the lifetime of a program.
    _realpaths = dict()

    def __init__(self, user_api=None, prefixes=None,  modules=None):
        if modules is None:
            self.prefixes = [] if prefixes is None else prefixes
            self.user_api = [] if user_api is None else user_api

            self.modules = []
            self._load_modules()
            self._warn_if_incompatible_openmp()
        else:
            self.modules = modules

    def get_modules(self, key, values):
        """Return all modules such that values contains module[key]"""
        if key == "user_api" and values is None:
            values = list(_ALL_USER_APIS)
        if not isinstance(values, list):
            values = [values]
        modules = [module for module in self.modules
                   if getattr(module, key) in values]
        return _ThreadpoolInfo(modules=modules)

    def todicts(self):
        """Return info as a list of dicts"""
        return [module.todict() for module in self.modules]

    def __len__(self):
        return len(self.modules)

    def __iter__(self):
        yield from self.modules

    def __eq__(self, other):
        return self.modules == other.modules

    def _load_modules(self):
        """Loop through loaded libraries and store supported ones"""
        if sys.platform == "darwin":
            self._find_modules_with_dyld()
        elif sys.platform == "win32":
            self._find_modules_with_enum_process_module_ex()
        else:
            self._find_modules_with_dl_iterate_phdr()

    def _find_modules_with_dl_iterate_phdr(self):
        """Loop through loaded libraries and return binders on supported ones

        This function is expected to work on POSIX system only.
        This code is adapted from code by Intel developper @anton-malakhov
        available at https://github.com/IntelPython/smp

        Copyright (c) 2017, Intel Corporation published under the BSD 3-Clause
        license
        """
        libc = self._get_libc()
        if not hasattr(libc, "dl_iterate_phdr"):  # pragma: no cover
            return []

        # Callback function for `dl_iterate_phdr` which is called for every
        # module loaded in the current process until it returns 1.
        def match_module_callback(info, size, data):
            # Get the path of the current module
            filepath = info.contents.dlpi_name
            if filepath:
                filepath = filepath.decode("utf-8")

                # Store the module if it is supported and selected
                self._make_module_from_path(filepath)
            return 0

        c_func_signature = ctypes.CFUNCTYPE(
            ctypes.c_int,  # Return type
            ctypes.POINTER(_dl_phdr_info), ctypes.c_size_t, ctypes.c_char_p)
        c_match_module_callback = c_func_signature(match_module_callback)

        data = ctypes.c_char_p(b"")
        libc.dl_iterate_phdr(c_match_module_callback, data)

    def _find_modules_with_dyld(self):
        """Loop through loaded libraries and return binders on supported ones

        This function is expected to work on OSX system only
        """
        libc = self._get_libc()
        if not hasattr(libc, "_dyld_image_count"):  # pragma: no cover
            return []

        n_dyld = libc._dyld_image_count()
        libc._dyld_get_image_name.restype = ctypes.c_char_p

        for i in range(n_dyld):
            filepath = ctypes.string_at(libc._dyld_get_image_name(i))
            filepath = filepath.decode("utf-8")

            # Store the module if it is supported and selected
            self._make_module_from_path(filepath)

    def _find_modules_with_enum_process_module_ex(self):
        """Loop through loaded libraries and return binders on supported ones

        This function is expected to work on windows system only.
        This code is adapted from code by Philipp Hagemeister @phihag available
        at https://stackoverflow.com/questions/17474574
        """
        from ctypes.wintypes import DWORD, HMODULE, MAX_PATH

        PROCESS_QUERY_INFORMATION = 0x0400
        PROCESS_VM_READ = 0x0010

        LIST_MODULES_ALL = 0x03

        ps_api = self._get_windll("Psapi")
        kernel_32 = self._get_windll("kernel32")

        h_process = kernel_32.OpenProcess(
            PROCESS_QUERY_INFORMATION | PROCESS_VM_READ,
            False, os.getpid())
        if not h_process:  # pragma: no cover
            raise OSError("Could not open PID %s" % os.getpid())

        try:
            buf_count = 256
            needed = DWORD()
            # Grow the buffer until it becomes large enough to hold all the
            # module headers
            while True:
                buf = (HMODULE * buf_count)()
                buf_size = ctypes.sizeof(buf)
                if not ps_api.EnumProcessModulesEx(
                        h_process, ctypes.byref(buf), buf_size,
                        ctypes.byref(needed), LIST_MODULES_ALL):
                    raise OSError("EnumProcessModulesEx failed")
                if buf_size >= needed.value:
                    break
                buf_count = needed.value // (buf_size // buf_count)

            count = needed.value // (buf_size // buf_count)
            h_modules = map(HMODULE, buf[:count])

            # Loop through all the module headers and get the module path
            buf = ctypes.create_unicode_buffer(MAX_PATH)
            n_size = DWORD()
            for h_module in h_modules:

                # Get the path of the current module
                if not ps_api.GetModuleFileNameExW(
                        h_process, h_module, ctypes.byref(buf),
                        ctypes.byref(n_size)):
                    raise OSError("GetModuleFileNameEx failed")
                filepath = buf.value

                # Store the module if it is supported and selected
                self._make_module_from_path(filepath)
        finally:
            kernel_32.CloseHandle(h_process)

    def _make_module_from_path(self, filepath):
        """Store a module if it is supported and selected"""
        # Required to resolve symlinks
        filepath = self._realpath(filepath)
        # `lower` required to take account of OpenMP dll case on Windows
        # (vcomp, VCOMP, Vcomp, ...)
        filename = os.path.basename(filepath).lower()

        # Loop through supported modules to find if this filename corresponds
        # to a supported module.
        for module_class, candidate_module in _SUPPORTED_MODULES.items():
            # check if filename matches a supported prefix
            prefix = self._check_prefix(filename,
                                        candidate_module["filename_prefixes"])

            # filename does not match any of the prefixes of the candidate
            # module. move to next module.
            if prefix is None:
                continue

            # filename matches a prefix. Check if it matches the request. If
            # so, create and store the module.
            user_api = candidate_module["user_api"]
            internal_api = candidate_module["internal_api"]
            if prefix in self.prefixes or user_api in self.user_api:
                module_class = globals()[module_class]
                module = module_class(filepath, prefix, user_api, internal_api)
                self.modules.append(module)

    def _check_prefix(self, library_basename, filename_prefixes):
        """Return the prefix library_basename starts with

        Return None if none matches.
        """
        for prefix in filename_prefixes:
            if library_basename.startswith(prefix):
                return prefix
        return None

    def _warn_if_incompatible_openmp(self):
        """Raise a warning if llvm-OpenMP and intel-OpenMP are both loaded"""
        if sys.platform != 'linux':
            # Only raise the warning on linux
            return

        prefixes = [module.prefix for module in self.modules]
        msg = textwrap.dedent(
            """
            Found Intel OpenMP ('libiomp') and LLVM OpenMP ('libomp') loaded at
            the same time. Both libraries are known to be incompatible and this
            can cause random crashes or deadlocks on Linux when loaded in the
            same Python program.
            Using threadpoolctl may cause crashes or deadlocks. For more
            information and possible workarounds, please see
                https://github.com/joblib/threadpoolctl/blob/master/multiple_openmp.md
            """)
        if 'libomp' in prefixes and 'libiomp' in prefixes:
            warnings.warn(msg, RuntimeWarning)

    @classmethod
    def _get_libc(cls):
        """Load the lib-C for unix systems."""
        libc = cls._system_libraries.get("libc")
        if libc is None:
            libc_name = find_library("c")
            if libc_name is None:  # pragma: no cover
                return None
            libc = ctypes.CDLL(libc_name, mode=_RTLD_NOLOAD)
            cls._system_libraries["libc"] = libc
        return libc

    @classmethod
    def _get_windll(cls, dll_name):
        """Load a windows DLL"""
        dll = cls._system_libraries.get(dll_name)
        if dll is None:
            dll = ctypes.WinDLL("{}.dll".format(dll_name))
            cls._system_libraries[dll_name] = dll
        return dll

    @classmethod
    def _realpath(cls, filepath, cache_limit=10000):
        """Small caching wrapper around os.path.realpath to limit system calls
        """
        rpath = cls._realpaths.get(filepath)
        if rpath is None:
            rpath = os.path.realpath(filepath)
            if len(cls._realpaths) < cache_limit:
                # If we drop support for Python 2.7, we could use
                # functools.lru_cache with maxsize=10000 instead.
                cls._realpaths[filepath] = rpath
        return rpath


@_format_docstring(
    USER_APIS=", ".join('"{}"'.format(api) for api in _ALL_USER_APIS),
    INTERNAL_APIS=", ".join('"{}"'.format(api) for api in _ALL_INTERNAL_APIS))
class _Module(ABC):
    """Abstract base class for the modules

    A module is represented by the following information:
      - "user_api" : user API. Possible values are {USER_APIS}.
      - "internal_api" : internal API. Possible values are {INTERNAL_APIS}.
      - "prefix" : prefix of the shared library's name.
      - "filepath" : path to the loaded module.
      - "version" : version of the library (if available).
      - "num_threads" : the current thread limit.

    In addition, each module may contain internal_api specific entries.
    """
    def __init__(self, filepath=None, prefix=None, user_api=None,
                 internal_api=None):
        self.filepath = filepath
        self.prefix = prefix
        self.user_api = user_api
        self.internal_api = internal_api
        self._dynlib = ctypes.CDLL(filepath, mode=_RTLD_NOLOAD)
        self.version = self.get_version()
        self.num_threads = self.get_num_threads()
        self._get_extra_info()

    def __eq__(self, other):
        return self.todict() == other.todict()

    def todict(self):
        """Return relevant info wrapped in a dict"""
        return {k: v for k, v in vars(self).items() if not k.startswith("_")}

    @abstractmethod
    def get_version(self):
        """Return the version of the shared library"""
        pass  # pragma: no cover

    @abstractmethod
    def get_num_threads(self):
        """Return the maximum number of threads available to use"""
        pass  # pragma: no cover

    @abstractmethod
    def set_num_threads(self, num_threads):
        """Set the maximum number of threads to use"""
        pass  # pragma: no cover

    @abstractmethod
    def _get_extra_info(self):
        """Add additional module specific information"""
        pass  # pragma: no cover


class _OpenBLASModule(_Module):
    """Module class for OpenBLAS"""
    def get_version(self):
        # None means OpenBLAS is not loaded or version < 0.3.4, since OpenBLAS
        # did not expose its version before that.
        get_config = getattr(self._dynlib, "openblas_get_config",
                             lambda: None)
        get_config.restype = ctypes.c_char_p
        config = get_config().split()
        if config[0] == b"OpenBLAS":
            return config[1].decode("utf-8")
        return None

    def get_num_threads(self):
        get_func = getattr(self._dynlib, "openblas_get_num_threads",
                           lambda: None)
        return get_func()

    def set_num_threads(self, num_threads):
        set_func = getattr(self._dynlib, "openblas_set_num_threads",
                           lambda num_threads: None)
        return set_func(num_threads)

    def _get_extra_info(self):
        self.threading_layer = self.get_threading_layer()

    def get_threading_layer(self):
        """Return the threading layer of OpenBLAS"""
        threading_layer = self._dynlib.openblas_get_parallel()
        if threading_layer == 2:
            return "openmp"
        elif threading_layer == 1:
            return "pthreads"
        return "disabled"


class _BLISModule(_Module):
    """Module class for BLIS"""
    def get_version(self):
        get_version_ = getattr(self._dynlib, "bli_info_get_version_str",
                               lambda: None)
        get_version_.restype = ctypes.c_char_p
        return get_version_().decode("utf-8")

    def get_num_threads(self):
        get_func = getattr(self._dynlib, "bli_thread_get_num_threads",
                           lambda: None)
        num_threads = get_func()
        # by default BLIS is single-threaded and get_num_threads
        # returns -1. We map it to 1 for consistency with other libraries.
        return 1 if num_threads == -1 else num_threads

    def set_num_threads(self, num_threads):
        set_func = getattr(self._dynlib, "bli_thread_set_num_threads",
                           lambda num_threads: None)
        return set_func(num_threads)

    def _get_extra_info(self):
        self.threading_layer = self.get_threading_layer()

    def get_threading_layer(self):
        """Return the threading layer of BLIS"""
        if self._dynlib.bli_info_get_enable_openmp():
            return "openmp"
        elif self._dynlib.bli_info_get_enable_pthreads():
            return "pthreads"
        return "disabled"


class _MKLModule(_Module):
    """Module class for MKL"""
    def get_version(self):
        res = ctypes.create_string_buffer(200)
        self._dynlib.mkl_get_version_string(res, 200)

        version = res.value.decode("utf-8")
        group = re.search(r"Version ([^ ]+) ", version)
        if group is not None:
            version = group.groups()[0]
        return version.strip()

    def get_num_threads(self):
        get_func = getattr(self._dynlib, "MKL_Get_Max_Threads", lambda: None)
        return get_func()

    def set_num_threads(self, num_threads):
        set_func = getattr(self._dynlib, "MKL_Set_Num_Threads",
                           lambda num_threads: None)
        return set_func(num_threads)

    def _get_extra_info(self):
        self.threading_layer = self.get_threading_layer()

    def get_threading_layer(self):
        """Return the threading layer of MKL"""
        # The function mkl_set_threading_layer returns the current threading
        # layer. Calling it with an invalid threading layer allows us to safely
        # get the threading layer
        set_threading_layer = getattr(self._dynlib, "MKL_Set_Threading_Layer",
                                      lambda layer: -1)
        layer_map = {0: "intel", 1: "sequential", 2: "pgi",
                     3: "gnu", 4: "tbb", -1: "not specified"}
        return layer_map[set_threading_layer(-1)]


class _OpenMPModule(_Module):
    """Module class for OpenMP"""
    def get_version(self):
        # There is no way to get the version number programmatically in OpenMP.
        return None

    def get_num_threads(self):
        get_func = getattr(self._dynlib, "omp_get_max_threads", lambda: None)
        return get_func()

    def set_num_threads(self, num_threads):
        set_func = getattr(self._dynlib, "omp_set_num_threads",
                           lambda num_threads: None)
        return set_func(num_threads)

    def _get_extra_info(self):
        pass
