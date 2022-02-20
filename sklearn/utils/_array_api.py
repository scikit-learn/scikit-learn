"""Tools to support array_api."""
import numpy
from .._config import get_config


class _ArrayAPIWrapper:
    def __init__(self, array_namespace):
        self._namespace = array_namespace

    def __getattr__(self, name):
        return getattr(self._namespace, name)

    def astype(self, x, dtype, *, copy=True, casting="unsafe"):
        # support casting for NumPy
        if self._namespace.__name__ == "numpy.array_api":
            x_np = numpy.asarray(x).astype(dtype, casting=casting, copy=copy)
            return self._namespace.asarray(x_np)

        f = self._namespace.astype
        return f(x, dtype, copy=copy)

    def asarray(self, obj, *, dtype=None, device=None, copy=None, order=None):
        f = self._namespace.asarray

        # support order in NumPy
        if self._namespace.__name__ == "numpy.array_api":
            if copy:
                x_np = numpy.array(obj, dtype=dtype, order=order, copy=True)
            else:
                x_np = numpy.asarray(obj, dtype=dtype, order=order)
            return f(x_np)
        return f(obj, dtype=dtype, device=device, copy=copy)

    def may_share_memory(self, a, b):
        # support may_share_memory in NumPy
        if self._namespace.__name__ == "numpy.array_api":
            return numpy.may_share_memory(a, b)

        # The safe choice is to return True for all other array_api Arrays
        return True


class _NumPyApiWrapper:
    def __getattr__(self, name):
        return getattr(numpy, name)

    def astype(self, x, dtype, *, copy=True, casting="unsafe"):
        # astype is not defined in the top level NumPy namespace
        return x.astype(dtype, copy=copy, casting=casting)

    def asarray(self, obj, *, dtype=None, device=None, copy=None, order=None):
        # copy is in the ArrayAPI spec but not in NumPy's asarray
        if copy:
            return numpy.array(obj, dtype=dtype, order=order, copy=True)
        else:
            return numpy.asarray(obj, dtype=dtype, order=order)

    def unique_inverse(self, x):
        return numpy.unique(x, return_inverse=True)

    def unique_counts(self, x):
        return numpy.unique(x, return_counts=True)

    def unique_values(self, x):
        return numpy.unique(x)

    def concat(self, arrays, *, axis=None):
        return numpy.concatenate(arrays, axis=axis)


def get_namespace(*arrays):
    """Get namespace of arrays.

    Parameters
    ----------
    *arrays : array objects
        Array objects.

    Returns
    -------
    namespace : module
        Namespace shared by array objects.
    """
    # `arrays` contains one or more arrays, or possibly Python scalars (accepting
    # those is a matter of taste, but doesn't seem unreasonable).
    # Returns a tuple: (array_namespace, is_array_api)

    if not get_config()["array_api_dispatch"]:
        return _NumPyApiWrapper(), False

    namespaces = {
        x.__array_namespace__() if hasattr(x, "__array_namespace__") else None
        for x in arrays
        if not isinstance(x, (bool, int, float, complex))
    }

    if not namespaces:
        # one could special-case np.ndarray above or use np.asarray here if
        # older numpy versions need to be supported.
        raise ValueError("Unrecognized array input")

    if len(namespaces) != 1:
        raise ValueError(f"Multiple namespaces for array inputs: {namespaces}")

    (xp,) = namespaces
    if xp is None:
        # Use numpy as default
        return _NumPyApiWrapper(), False

    return _ArrayAPIWrapper(xp), True
