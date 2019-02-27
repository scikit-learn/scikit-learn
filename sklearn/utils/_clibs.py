"""
This module provides utilities to load C-libraries that relies on thread
pools and limit the maximal number of thread that can be used.
"""

# This code is adapted from code by Thomas Moreau <tomMoral> available at
# https://github.com/tomMoral/loky


import sys
import os
import threading
import ctypes
from ctypes.util import find_library
from contextlib import contextmanager


# Structure to cast the info on dynamically loaded library. See
# https://linux.die.net/man/3/dl_iterate_phdr for more details.
UINT_SYSTEM = ctypes.c_uint64 if sys.maxsize > 2**32 else ctypes.c_uint32
UINT_HALF_SYSTEM = ctypes.c_uint32 if sys.maxsize > 2**32 else ctypes.c_uint16


class dl_phdr_info(ctypes.Structure):
    _fields_ = [
        ("dlpi_addr",  UINT_SYSTEM),       # Base address of object
        ("dlpi_name",  ctypes.c_char_p),   # path to the library
        ("dlpi_phdr",  ctypes.c_void_p),   # pointer on dlpi_headers
        ("dlpi_phnum",  UINT_HALF_SYSTEM)  # number of element in dlpi_phdr
        ]


class _CLibsWrapper:
    # Wrapper around classic C-libraries for scientific computations to set and
    # get the maximum number of threads they are allowed to used for inner
    # parallelism.

    # Supported C-libraries for this wrapper, index with their name. The items
    # hold the name of the library file and the functions to call.
    SUPPORTED_CLIBS = {
        "openmp_intel": (
            "libiomp", "omp_set_num_threads", "omp_get_max_threads"),
        "openmp_gnu": (
            "libgomp", "omp_set_num_threads", "omp_get_max_threads"),
        "openmp_llvm": (
            "libomp", "omp_set_num_threads", "omp_get_max_threads"),
        "openmp_win32": (
            "vcomp", "omp_set_num_threads", "omp_get_max_threads"),
        "openblas": (
            "libopenblas", "openblas_set_num_threads",
            "openblas_get_num_threads"),
        "mkl": (
            "libmkl_rt", "MKL_Set_Num_Threads", "MKL_Get_Max_Threads"),
        "mkl_win32": (
            "mkl_rt", "MKL_Set_Num_Threads", "MKL_Get_Max_Threads")}

    cls_thread_locals = threading.local()

    def __init__(self):
        self._load()

    def _load(self):
        for clib, (module_name, _, _) in self.SUPPORTED_CLIBS.items():
            setattr(self, clib, self._load_lib(module_name))

    def _unload(self):
        for clib, (module_name, _, _) in self.SUPPORTED_CLIBS.items():
            delattr(self, clib)

    def set_thread_limits(self, limits=1, subset=None):
        """Limit maximal number of threads used by supported C-libraries"""
        if isinstance(limits, int):
            if subset in ("all", None):
                clibs = self.SUPPORTED_CLIBS.keys()
            elif subset == "blas":
                clibs = ("openblas", "mkl", "mkl_win32")
            elif subset == "openmp":
                clibs = (c for c in self.SUPPORTED_CLIBS if "openmp" in c)
            else:
                raise ValueError("subset must be either 'all', 'blas' or "
                                 "'openmp'. Got {} instead.".format(subset))
            limits = {clib: limits for clib in clibs}

        if not isinstance(limits, dict):
            raise TypeError("limits must either be an int or a dict. Got {} "
                            "instead".format(type(limits)))

        dynamic_threadpool_size = {}
        self._load()
        for clib, (_, _set, _) in self.SUPPORTED_CLIBS.items():
            if clib in limits:
                module = getattr(self, clib, None)
                if module is not None:
                    _set = getattr(module, _set)
                    _set(limits[clib])
                    dynamic_threadpool_size[clib] = True
                else:
                    dynamic_threadpool_size[clib] = False
            else:
                dynamic_threadpool_size[clib] = False
        self._unload()
        return dynamic_threadpool_size

    def get_thread_limits(self):
        """Return maximal number of threads available for supported C-libraries
        """
        limits = {}
        self._load()
        for clib, (_, _, _get) in self.SUPPORTED_CLIBS.items():
            module = getattr(self, clib, None)
            if module is not None:
                _get = getattr(module, _get)
                limits[clib] = _get()
            else:
                limits[clib] = None
        self._unload()
        return limits

    def get_openblas_version(self):
        module = getattr(self, "openblas", None)
        if module is not None:
            get_config = getattr(module, "openblas_get_config")
            get_config.restype = ctypes.c_char_p
            config = get_config().split()
            if config[0] == b"OpenBLAS":
                return config[1].decode('utf-8')
            return
        return

    def _load_lib(self, module_name):
        """Return a binder on module_name by looping through loaded libraries
        """
        if sys.platform == "darwin":
            return self._find_with_clibs_dyld(module_name)
        elif sys.platform == "win32":
            return self._find_with_clibs_enum_process_module_ex(module_name)
        return self._find_with_clibs_dl_iterate_phdr(module_name)

    def _find_with_clibs_dl_iterate_phdr(self, module_name):
        """Return a binder on module_name by looping through loaded libraries

        This function is expected to work on POSIX system only.
        This code is adapted from code by Intel developper @anton-malakhov
        available at https://github.com/IntelPython/smp

        Copyright (c) 2017, Intel Corporation published under the BSD 3-Clause
        license
        """
        self.cls_thread_locals._module_path = None

        libc = self._get_libc()
        if not hasattr(libc, "dl_iterate_phdr"):
            return

        # Callback function for `dl_iterate_phdr` which is called for every
        # module loaded in the current process until it returns 1.
        def match_module_callback(info, size, module_name):

            # recast the name of the module as a string
            module_name = ctypes.string_at(module_name).decode('utf-8')

            # Get the name of the current library
            module_path = info.contents.dlpi_name

            # If the current library is the one we are looking for, store the
            # path and return 1 to stop the loop in `dl_iterate_phdr`.
            if module_path:
                module_path = module_path.decode("utf-8")
                if os.path.basename(module_path).startswith(module_name):
                    self.cls_thread_locals._module_path = module_path
                    return 1
            return 0

        c_func_signature = ctypes.CFUNCTYPE(
            ctypes.c_int,  # Return type
            ctypes.POINTER(dl_phdr_info), ctypes.c_size_t, ctypes.c_char_p)
        c_match_module_callback = c_func_signature(match_module_callback)

        data = ctypes.c_char_p(module_name.encode('utf-8'))
        res = libc.dl_iterate_phdr(c_match_module_callback, data)
        if res == 1:
            return ctypes.CDLL(self.cls_thread_locals._module_path)

    def _find_with_clibs_dyld(self, module_name):
        """Return a binder on module_name by looping through loaded libraries

        This function is expected to work on OSX system only
        """
        libc = self._get_libc()
        if not hasattr(libc, "_dyld_image_count"):
            return

        found_module_path = None

        n_dyld = libc._dyld_image_count()
        libc._dyld_get_image_name.restype = ctypes.c_char_p

        for i in range(n_dyld):
            module_path = ctypes.string_at(libc._dyld_get_image_name(i))
            module_path = module_path.decode("utf-8")
            if os.path.basename(module_path).startswith(module_name):
                found_module_path = module_path

        if found_module_path:
            return ctypes.CDLL(found_module_path)

    def _find_with_clibs_enum_process_module_ex(self, module_name):
        """Return a binder on module_name by looping through loaded libraries

        This function is expected to work on windows system only.
        This code is adapted from code by Philipp Hagemeister @phihag available
        at https://stackoverflow.com/questions/17474574
        """
        from ctypes.wintypes import DWORD, HMODULE, MAX_PATH

        PROCESS_QUERY_INFORMATION = 0x0400
        PROCESS_VM_READ = 0x0010

        LIST_MODULES_ALL = 0x03

        Psapi = self._get_windll('Psapi')
        Kernel32 = self._get_windll('kernel32')

        hProcess = Kernel32.OpenProcess(
            PROCESS_QUERY_INFORMATION | PROCESS_VM_READ,
            False, os.getpid())
        if not hProcess:
            raise OSError('Could not open PID %s' % os.getpid())

        found_module_path = None
        try:
            buf_count = 256
            needed = DWORD()
            # Grow the buffer until it becomes large enough to hold all the
            # module headers
            while True:
                buf = (HMODULE * buf_count)()
                buf_size = ctypes.sizeof(buf)
                if not Psapi.EnumProcessModulesEx(
                        hProcess, ctypes.byref(buf), buf_size,
                        ctypes.byref(needed), LIST_MODULES_ALL):
                    raise OSError('EnumProcessModulesEx failed')
                if buf_size >= needed.value:
                    break
                buf_count = needed.value // (buf_size // buf_count)

            count = needed.value // (buf_size // buf_count)
            hModules = map(HMODULE, buf[:count])

            # Loop through all the module headers and get the module file name
            buf = ctypes.create_unicode_buffer(MAX_PATH)
            nSize = DWORD()
            for hModule in hModules:
                if not Psapi.GetModuleFileNameExW(
                        hProcess, hModule, ctypes.byref(buf),
                        ctypes.byref(nSize)):
                    raise OSError('GetModuleFileNameEx failed')
                module_path = buf.value
                module_basename = os.path.basename(module_path).lower()
                if module_basename.startswith(module_name):
                    found_module_path = module_path
        finally:
            Kernel32.CloseHandle(hProcess)

        if found_module_path:
            return ctypes.CDLL(found_module_path)

    def _get_libc(self):
        if not hasattr(self, "libc"):
            libc_name = find_library("c")
            if libc_name is None:
                self.libc = None
            self.libc = ctypes.CDLL(libc_name)

        return self.libc

    def _get_windll(self, dll_name):
        if not hasattr(self, dll_name):
            setattr(self, dll_name, ctypes.WinDLL("{}.dll".format(dll_name)))

        return getattr(self, dll_name)


_clibs_wrapper = None


def _get_wrapper(reload_clib=False):
    """Helper function to only create one wrapper per thread."""
    global _clibs_wrapper
    if _clibs_wrapper is None:
        _clibs_wrapper = _CLibsWrapper()
    if reload_clib:
        _clibs_wrapper._load()
    return _clibs_wrapper


def set_thread_limits(limits=1, subset=None, reload_clib=False):
    """Limit the number of threads available for threadpools in supported C-lib

    Set the maximal number of thread that can be used in thread pools used in
    the supported C-libraries. This function works for libraries that are
    already loaded in the interpreter and can be changed dynamically.

    Parameters
    ----------
    limits : int or dict, (default=1)
        Maximum number of thread that can be used in thread pools

        If int, sets the maximum number of thread to `limits` for each C-lib
        selected by `subset`.

        If dict(supported_libraries: max_threads), sets a custom maximum number
        of thread for each C-lib.

    subset : string or None, optional (default="all")
        Subset of C-libs to limit. Used only if `limits` is an int

        "all" : limit all supported C-libs.

        "blas" : limit only BLAS supported C-libs.

        "openmp" : limit only OpenMP supported C-libs. It can affect the number
                   of threads used by the BLAS C-libs if they rely on OpenMP.

    reload_clib : bool, (default=False)
        If `reload_clib` is `True`, first loop through the loaded libraries to
        ensure that this function is called on all available libraries.

    Returns
    -------
    dynamic_threadpool_size : dict
        contains pairs `('clib': boolean)` which are True if `clib` have been
        found and can be used to scale the maximal number of threads
        dynamically.
    """
    wrapper = _get_wrapper(reload_clib)
    return wrapper.set_thread_limits(limits, subset)


def get_thread_limits(reload_clib=True):
    """Return maximal thread number for threadpools in supported C-lib

    Parameters
    ----------
    reload_clib : bool, (default=True)
        If `reload_clib` is `True`, first loop through the loaded libraries to
        ensure that this function is called on all available libraries.

    Returns
    -------
    thread_limits : dict
        Contains the maximal number of threads that can be used in supported
        libraries or None when the library is not available. The key of the
        dictionary are "openmp_gnu", "openmp_intel", "openmp_win32",
        "openmp_llvm", "openblas", "mkl" and "mkl_win32".
    """
    wrapper = _get_wrapper(reload_clib)
    return wrapper.get_thread_limits()


@contextmanager
def thread_limits_context(limits=1, subset=None):
    """Context manager for C-libs thread limits

    Parameters
    ----------
    limits : int or dict, (default=1)
        Maximum number of thread that can be used in thread pools

        If int, sets the maximum number of thread to `limits` for each C-lib
        selected by `subset`.

        If dict(supported_libraries: max_threads), sets a custom maximum number
        of thread for each C-lib.

    subset : string or None, optional (default="all")
        Subset of C-libs to limit. Used only if `limits` is an int

        "all" : limit all supported C-libs.

        "blas" : limit only BLAS supported C-libs.

        "openmp" : limit only OpenMP supported C-libs. It can affect the number
                   of threads used by the BLAS C-libs if they rely on OpenMP.
    """
    old_limits = get_thread_limits()
    set_thread_limits(limits=limits, subset=subset)

    try:
        yield
    finally:
        set_thread_limits(limits=old_limits)


def get_openblas_version(reload_clib=True):
    """Return the OpenBLAS version

    Parameters
    ----------
    reload_clib : bool, (default=True)
        If `reload_clib` is `True`, first loop through the loaded libraries to
        ensure that this function is called on all available libraries.

    Returns
    -------
    version : string or None
        None means OpenBLAS is not loaded or version < 0.3.4, since OpenBLAS
        did not expose it's verion before that.
    """
    wrapper = _get_wrapper(reload_clib)
    return wrapper.get_openblas_version()
