__all__ = [
    'load_ctypes_library',
    'addressof_array',
    'array_as_ctype'
    ]

import numpy as N
import ctypes

def load_ctypes_library(libname, loader_path):
    import sys, os
    if sys.platform == 'win32':
        libname = '%s.dll' % libname
    else:
        libname = '%s.so' % libname
    loader_path = os.path.abspath(loader_path)
    if not os.path.isdir(loader_path):
        libdir = os.path.dirname(loader_path)
    else:
        libdir = loader_path
    libpath = os.path.join(libdir, libname)
    if sys.platform == 'win32':
        return ctypes.cdll.load(libpath)
    else:
        return ctypes.cdll.LoadLibrary(libpath)

def addressof_array(a):
    return N.cast[N.intp](int(a.__array_data__[0], 16))

def array_as_ctype(a, ctype):
    return ctypes.cast(addressof_array(a), ctypes.POINTER(ctype))
