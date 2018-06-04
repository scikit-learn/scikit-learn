from __future__ import absolute_import, division, print_function

import platform
import sys
import types
import warnings


PY2 = sys.version_info[0] == 2
PYPY = platform.python_implementation() == "PyPy"


if PYPY or sys.version_info[:2] >= (3, 6):
    ordered_dict = dict
else:
    from collections import OrderedDict
    ordered_dict = OrderedDict


if PY2:
    from UserDict import IterableUserDict

    # We 'bundle' isclass instead of using inspect as importing inspect is
    # fairly expensive (order of 10-15 ms for a modern machine in 2016)
    def isclass(klass):
        return isinstance(klass, (type, types.ClassType))

    # TYPE is used in exceptions, repr(int) is different on Python 2 and 3.
    TYPE = "type"

    def iteritems(d):
        return d.iteritems()

    # Python 2 is bereft of a read-only dict proxy, so we make one!
    class ReadOnlyDict(IterableUserDict):
        """
        Best-effort read-only dict wrapper.
        """

        def __setitem__(self, key, val):
            # We gently pretend we're a Python 3 mappingproxy.
            raise TypeError("'mappingproxy' object does not support item "
                            "assignment")

        def update(self, _):
            # We gently pretend we're a Python 3 mappingproxy.
            raise AttributeError("'mappingproxy' object has no attribute "
                                 "'update'")

        def __delitem__(self, _):
            # We gently pretend we're a Python 3 mappingproxy.
            raise TypeError("'mappingproxy' object does not support item "
                            "deletion")

        def clear(self):
            # We gently pretend we're a Python 3 mappingproxy.
            raise AttributeError("'mappingproxy' object has no attribute "
                                 "'clear'")

        def pop(self, key, default=None):
            # We gently pretend we're a Python 3 mappingproxy.
            raise AttributeError("'mappingproxy' object has no attribute "
                                 "'pop'")

        def popitem(self):
            # We gently pretend we're a Python 3 mappingproxy.
            raise AttributeError("'mappingproxy' object has no attribute "
                                 "'popitem'")

        def setdefault(self, key, default=None):
            # We gently pretend we're a Python 3 mappingproxy.
            raise AttributeError("'mappingproxy' object has no attribute "
                                 "'setdefault'")

        def __repr__(self):
            # Override to be identical to the Python 3 version.
            return "mappingproxy(" + repr(self.data) + ")"

    def metadata_proxy(d):
        res = ReadOnlyDict()
        res.data.update(d)  # We blocked update, so we have to do it like this.
        return res

else:
    def isclass(klass):
        return isinstance(klass, type)

    TYPE = "class"

    def iteritems(d):
        return d.items()

    def metadata_proxy(d):
        return types.MappingProxyType(dict(d))


def import_ctypes():
    """
    Moved into a function for testability.
    """
    import ctypes
    return ctypes


if not PY2:
    def just_warn(*args, **kw):
        """
        We only warn on Python 3 because we are not aware of any concrete
        consequences of not setting the cell on Python 2.
        """
        warnings.warn(
            "Missing ctypes.  Some features like bare super() or accessing "
            "__class__ will not work with slots classes.",
            RuntimeWarning,
            stacklevel=2,
        )
else:
    def just_warn(*args, **kw):  # pragma: nocover
        """
        We only warn on Python 3 because we are not aware of any concrete
        consequences of not setting the cell on Python 2.
        """


def make_set_closure_cell():
    """
    Moved into a function for testability.
    """
    if PYPY:  # pragma: no cover
        def set_closure_cell(cell, value):
            cell.__setstate__((value,))
    else:
        try:
            ctypes = import_ctypes()

            set_closure_cell = ctypes.pythonapi.PyCell_Set
            set_closure_cell.argtypes = (ctypes.py_object, ctypes.py_object)
            set_closure_cell.restype = ctypes.c_int
        except Exception:
            # We try best effort to set the cell, but sometimes it's not
            # possible.  For example on Jython or on GAE.
            set_closure_cell = just_warn
    return set_closure_cell


set_closure_cell = make_set_closure_cell()
