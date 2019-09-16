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
import warnings
from ctypes.util import find_library

__version__ = '1.1.0'
__all__ = ["threadpool_limits", "threadpool_info"]

# Cache for libc under POSIX and a few system libraries under Windows
_system_libraries = {}

# Cache for calls to os.path.realpath on system libraries to reduce the
# impact of slow system calls (e.g. stat) on slow filesystem
_realpaths = dict()

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
        ("dlpi_name",  ctypes.c_char_p),   # path to the library
        ("dlpi_phdr",  ctypes.c_void_p),   # pointer on dlpi_headers
        ("dlpi_phnum",  _SYSTEM_UINT_HALF)  # number of element in dlpi_phdr
    ]


# List of the supported implementations. The items hold the prefix of loaded
# shared objects, the name of the internal_api to call, matching the
# MAP_API_TO_FUNC keys and the name of the user_api, in {"blas", "openmp"}.

_SUPPORTED_IMPLEMENTATIONS = [
    {
        "user_api": "openmp",
        "internal_api": "openmp",
        "filename_prefixes": ("libiomp", "libgomp", "libomp", "vcomp",),
    },
    {
        "user_api": "blas",
        "internal_api": "openblas",
        "filename_prefixes": ("libopenblas",),
    },
    {
        "user_api": "blas",
        "internal_api": "mkl",
        "filename_prefixes": ("libmkl_rt", "mkl_rt",),
    },
    {
        "user_api": "blas",
        "internal_api": "blis",
        "filename_prefixes": ("libblis",),
    },
]

# map a internal_api (openmp, openblas, mkl) to set and get functions

_MAP_API_TO_FUNC = {
    "openmp": {
        "set_num_threads": "omp_set_num_threads",
        "get_num_threads": "omp_get_max_threads"},
    "openblas": {
        "set_num_threads": "openblas_set_num_threads",
        "get_num_threads": "openblas_get_num_threads"},
    "mkl": {
        "set_num_threads": "MKL_Set_Num_Threads",
        "get_num_threads": "MKL_Get_Max_Threads"},
    "blis": {
        "set_num_threads": "bli_thread_set_num_threads",
        "get_num_threads": "bli_thread_get_num_threads"}
}

# Helpers for the doc and test names

_ALL_USER_APIS = set(impl['user_api'] for impl in _SUPPORTED_IMPLEMENTATIONS)
_ALL_PREFIXES = [prefix
                 for impl in _SUPPORTED_IMPLEMENTATIONS
                 for prefix in impl['filename_prefixes']]
_ALL_INTERNAL_APIS = list(_MAP_API_TO_FUNC.keys())


def _realpath(filepath, cache_limit=10000):
    """Small caching wrapper around os.path.realpath to limit system calls"""
    rpath = _realpaths.get(filepath)
    if rpath is None:
        rpath = os.path.realpath(filepath)
        if len(_realpaths) < cache_limit:
            # If we drop support for Python 2.7, we could use functools.lru_cache
            # with maxsize=10000 instead.
            _realpaths[filepath] = rpath
    return rpath


def _format_docstring(*args, **kwargs):
    def decorator(o):
        o.__doc__ = o.__doc__.format(*args, **kwargs)
        return o

    return decorator


def _get_limit(prefix, user_api, limits):
    if prefix in limits:
        return limits[prefix]
    else:
        return limits[user_api]


@_format_docstring(ALL_PREFIXES=_ALL_PREFIXES,
                   INTERNAL_APIS=_ALL_INTERNAL_APIS)
def _set_threadpool_limits(limits, user_api=None):
    """Limit the maximal number of threads for threadpools in supported libs

    Set the maximal number of threads that can be used in thread pools used in
    the supported native libraries to `limit`. This function works for
    libraries that are already loaded in the interpreter and can be changed
    dynamically.

    The `limits` parameter can be either an integer or a dict to specify the
    maximal number of thread that can be used in thread pools. If it is an
    integer, sets the maximum number of thread to `limits` for each library
    selected by `user_api`. If it is a dictionary `{{key: max_threads}}`, this
    function sets a custom maximum number of thread for each `key` which can be
    either a `user_api` or a `prefix` for a specific library.

    The `user_api` parameter selects particular APIs of libraries to limit.
    Used only if `limits` is an int. If it is None, this function will apply to
    all supported libraries. If it is "blas", it will limit only BLAS supported
    libraries and if it is "openmp", only OpenMP supported libraries will be
    limited. Note that the latter can affect the number of threads used by the
    BLAS libraries if they rely on OpenMP.

    Return a list with all the supported modules that have been found. Each
    module is represented by a dict with the following information:
      - 'filename_prefixes' : possible prefixes for the given internal_api.
            Possible values are {ALL_PREFIXES}.
      - 'prefix' : prefix of the specific implementation of this module.
      - 'internal_api': internal API.s Possible values are {INTERNAL_APIS}.
      - 'filepath': path to the loaded module.
      - 'version': version of the library implemented (if available).
      - 'num_threads': the theadpool size limit before changing it.
      - 'set_num_threads': callable to set the maximum number of threads
      - 'get_num_threads': callable to get the current number of threads
      - 'dynlib': the instance of ctypes.CDLL use to access the dynamic
        library.
    """
    if isinstance(limits, int):
        if user_api is None:
            user_api = _ALL_USER_APIS
        elif user_api in _ALL_USER_APIS:
            user_api = (user_api,)
        else:
            raise ValueError("user_api must be either in {} or None. Got {} "
                             "instead.".format(_ALL_USER_APIS, user_api))
        limits = {api: limits for api in user_api}
        prefixes = []
    else:
        if isinstance(limits, list):
            # This should be a list of module, for compatibility with
            # the result from threadpool_info.
            limits = {module['prefix']: module['num_threads']
                      for module in limits}

        if not isinstance(limits, dict):
            raise TypeError("limits must either be an int, a list or a dict."
                            " Got {} instead".format(type(limits)))

        # With a dictionary, can set both specific limit for given modules
        # and global limit for user_api. Fetch each separately.
        prefixes = [module for module in limits if module in _ALL_PREFIXES]
        user_api = [module for module in limits if module in _ALL_USER_APIS]

    modules = _load_modules(prefixes=prefixes, user_api=user_api)
    for module in modules:
        # Workaround clang bug (TODO: report it)
        module['get_num_threads']()

    for module in modules:
        module['num_threads'] = module['get_num_threads']()
        num_threads = _get_limit(module['prefix'], module['user_api'], limits)
        if num_threads is not None:
            set_func = module['set_num_threads']
            set_func(num_threads)

    return modules


@_format_docstring(INTERNAL_APIS=_ALL_INTERNAL_APIS)
def threadpool_info():
    """Return the maximal number of threads for each detected library.

    Return a list with all the supported modules that have been found. Each
    module is represented by a dict with the following information:
      - 'prefix' : filename prefix of the specific implementation.
      - 'filepath': path to the loaded module.
      - 'internal_api': internal API. Possible values are {INTERNAL_APIS}.
      - 'version': version of the library implemented (if available).
      - 'num_threads': the current thread limit.
    """
    infos = []
    modules = _load_modules(user_api=_ALL_USER_APIS)
    for module in modules:
        module['num_threads'] = module['get_num_threads']()
        # by default BLIS is single-threaded and get_num_threads returns -1.
        # we map it to 1 for consistency with other libraries.
        if module['num_threads'] == -1 and module['internal_api'] == 'blis':
            module['num_threads'] = 1
        # Remove the wrapper for the module and its function
        del module['set_num_threads'], module['get_num_threads']
        del module['dynlib']
        del module['filename_prefixes']
        infos.append(module)
    return infos


def _get_version(dynlib, internal_api):
    if internal_api == "mkl":
        return _get_mkl_version(dynlib)
    elif internal_api == "openmp":
        # There is no way to get the version number programmatically in
        # OpenMP.
        return None
    elif internal_api == "openblas":
        return _get_openblas_version(dynlib)
    elif internal_api == "blis":
        return _get_blis_version(dynlib)
    else:
        raise NotImplementedError("Unsupported API {}".format(internal_api))


def _get_mkl_version(mkl_dynlib):
    """Return the MKL version"""
    res = ctypes.create_string_buffer(200)
    mkl_dynlib.mkl_get_version_string(res, 200)

    version = res.value.decode('utf-8')
    group = re.search(r"Version ([^ ]+) ", version)
    if group is not None:
        version = group.groups()[0]
    return version.strip()


def _get_openblas_version(openblas_dynlib):
    """Return the OpenBLAS version

    None means OpenBLAS is not loaded or version < 0.3.4, since OpenBLAS
    did not expose its version before that.
    """
    get_config = getattr(openblas_dynlib, "openblas_get_config")
    get_config.restype = ctypes.c_char_p
    config = get_config().split()
    if config[0] == b"OpenBLAS":
        return config[1].decode('utf-8')
    return None


def _get_blis_version(blis_dynlib):
    """Return the BLIS version"""
    get_version = getattr(blis_dynlib, "bli_info_get_version_str")
    get_version.restype = ctypes.c_char_p
    return get_version().decode('utf-8')


# Loading utilities for dynamically linked shared objects

def _load_modules(prefixes=None, user_api=None):
    """Loop through loaded libraries and return supported ones."""
    if prefixes is None:
        prefixes = []
    if user_api is None:
        user_api = []
    if sys.platform == "darwin":
        return _find_modules_with_dyld(prefixes=prefixes, user_api=user_api)
    elif sys.platform == "win32":
        return _find_modules_with_enum_process_module_ex(
            prefixes=prefixes, user_api=user_api)
    else:
        return _find_modules_with_dl_iterate_phdr(
            prefixes=prefixes, user_api=user_api)


def _check_prefix(library_basename, filename_prefixes):
    """Return the prefix library_basename starts with or None if none matches
    """
    for prefix in filename_prefixes:
        if library_basename.startswith(prefix):
            return prefix
    return None


def _match_module(module_info, prefix, prefixes, user_api):
    """Return True if this module should be selected."""
    return prefix is not None and (prefix in prefixes or
                                   module_info['user_api'] in user_api)


def _make_module_info(filepath, module_info, prefix):
    """Make a dict with the information from the module."""
    filepath = os.path.normpath(filepath)
    dynlib = ctypes.CDLL(filepath)
    internal_api = module_info['internal_api']
    set_func = getattr(dynlib,
                       _MAP_API_TO_FUNC[internal_api]['set_num_threads'],
                       lambda num_threads: None)
    get_func = getattr(dynlib,
                       _MAP_API_TO_FUNC[internal_api]['get_num_threads'],
                       lambda: None)
    module_info = module_info.copy()
    module_info.update(dynlib=dynlib, filepath=filepath, prefix=prefix,
                       set_num_threads=set_func, get_num_threads=get_func,
                       version=_get_version(dynlib, internal_api))
    return module_info


def _get_module_info_from_path(filepath, prefixes, user_api, modules):
    # Required to resolve symlinks
    filepath =_realpath(filepath)
    # `lower` required to take account of OpenMP dll case on Windows
    # (vcomp, VCOMP, Vcomp, ...)
    filename = os.path.basename(filepath).lower()
    for info in _SUPPORTED_IMPLEMENTATIONS:
        prefix = _check_prefix(filename, info['filename_prefixes'])
        if _match_module(info, prefix, prefixes, user_api):
            modules.append(_make_module_info(filepath, info, prefix))


def _find_modules_with_dl_iterate_phdr(prefixes, user_api):
    """Loop through loaded libraries and return binders on supported ones

    This function is expected to work on POSIX system only.
    This code is adapted from code by Intel developper @anton-malakhov
    available at https://github.com/IntelPython/smp

    Copyright (c) 2017, Intel Corporation published under the BSD 3-Clause
    license
    """
    libc = _get_libc()
    if not hasattr(libc, "dl_iterate_phdr"):  # pragma: no cover
        return []

    _modules = []

    # Callback function for `dl_iterate_phdr` which is called for every
    # module loaded in the current process until it returns 1.
    def match_module_callback(info, size, data):
        # Get the path of the current module
        filepath = info.contents.dlpi_name
        if filepath:
            filepath = filepath.decode("utf-8")

            # Store the module in cls_thread_locals._module if it is
            # supported and selected
            _get_module_info_from_path(filepath, prefixes, user_api,
                                       _modules)
        return 0

    c_func_signature = ctypes.CFUNCTYPE(
        ctypes.c_int,  # Return type
        ctypes.POINTER(_dl_phdr_info), ctypes.c_size_t, ctypes.c_char_p)
    c_match_module_callback = c_func_signature(match_module_callback)

    data = ctypes.c_char_p(b'')
    libc.dl_iterate_phdr(c_match_module_callback, data)

    return _modules


def _find_modules_with_dyld(prefixes, user_api):
    """Loop through loaded libraries and return binders on supported ones

    This function is expected to work on OSX system only
    """
    libc = _get_libc()
    if not hasattr(libc, "_dyld_image_count"):  # pragma: no cover
        return []

    _modules = []

    n_dyld = libc._dyld_image_count()
    libc._dyld_get_image_name.restype = ctypes.c_char_p

    for i in range(n_dyld):
        filepath = ctypes.string_at(libc._dyld_get_image_name(i))
        filepath = filepath.decode("utf-8")

        # Store the module in cls_thread_locals._module if it is supported and
        # selected
        _get_module_info_from_path(filepath, prefixes, user_api, _modules)

    return _modules


def _find_modules_with_enum_process_module_ex(prefixes, user_api):
    """Loop through loaded libraries and return binders on supported ones

    This function is expected to work on windows system only.
    This code is adapted from code by Philipp Hagemeister @phihag available
    at https://stackoverflow.com/questions/17474574
    """
    from ctypes.wintypes import DWORD, HMODULE, MAX_PATH

    PROCESS_QUERY_INFORMATION = 0x0400
    PROCESS_VM_READ = 0x0010

    LIST_MODULES_ALL = 0x03

    ps_api = _get_windll('Psapi')
    kernel_32 = _get_windll('kernel32')

    h_process = kernel_32.OpenProcess(
        PROCESS_QUERY_INFORMATION | PROCESS_VM_READ,
        False, os.getpid())
    if not h_process:  # pragma: no cover
        raise OSError('Could not open PID %s' % os.getpid())

    _modules = []
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
                raise OSError('EnumProcessModulesEx failed')
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
                raise OSError('GetModuleFileNameEx failed')
            filepath = buf.value

            # Store the module in cls_thread_locals._module if it is
            # supported and selected
            _get_module_info_from_path(filepath, prefixes, user_api,
                                       _modules)
    finally:
        kernel_32.CloseHandle(h_process)

    return _modules


def _get_libc():
    """Load the lib-C for unix systems."""
    libc = _system_libraries.get("libc")
    if libc is None:
        libc_name = find_library("c")
        if libc_name is None:  # pragma: no cover
            return None
        libc = ctypes.CDLL(libc_name)
        _system_libraries["libc"] = libc
    return libc


def _get_windll(dll_name):
    """Load a windows DLL"""
    dll = _system_libraries.get(dll_name)
    if dll is None:
        dll = ctypes.WinDLL("{}.dll".format(dll_name))
        _system_libraries[dll_name] = dll
    return dll


class threadpool_limits:
    """Change the maximal number of threads that can be used in thread pools.

    This class can be used either as a function (the construction of this
    object limits the number of threads) or as a context manager, in a `with`
    block.

    Set the maximal number of threads that can be used in thread pools used in
    the supported libraries to `limit`. This function works for libraries that
    are already loaded in the interpreter and can be changed dynamically.

    The `limits` parameter can be either an integer or a dict to specify the
    maximal number of thread that can be used in thread pools. If it is an
    integer, sets the maximum number of thread to `limits` for each library
    selected by `user_api`. If it is a dictionary `{{key: max_threads}}`, this
    function sets a custom maximum number of thread for each `key` which can be
    either a `user_api` or a `prefix` for a specific library. If None, this
    function does not do anything.

    The `user_api` parameter selects particular APIs of libraries to limit.
    Used only if `limits` is an int. If it is None, this function will apply to
    all supported libraries. If it is "blas", it will limit only BLAS supported
    libraries and if it is "openmp", only OpenMP supported libraries will be
    limited. Note that the latter can affect the number of threads used by the
    BLAS libraries if they rely on OpenMP.
    """
    def __init__(self, limits=None, user_api=None):
        self._user_api = _ALL_USER_APIS if user_api is None else [user_api]

        if limits is not None:
            self._original_limits = _set_threadpool_limits(
                limits=limits, user_api=user_api)
        else:
            self._original_limits = None

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.unregister()

    def unregister(self):
        if self._original_limits is not None:
            for module in self._original_limits:
                module['set_num_threads'](module['num_threads'])

    def get_original_num_threads(self):
        original_limits = self._original_limits or threadpool_info()

        num_threads = {}
        warning_apis = []

        for user_api in self._user_api:
            limits = [module['num_threads'] for module in original_limits
                      if module['user_api'] == user_api]
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
            warnings.warn("Multiple value possible for following user apis: "
                          + ', '.join(warning_apis) + ". Returning the minimum.")

        return num_threads
