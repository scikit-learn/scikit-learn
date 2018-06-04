from __future__ import absolute_import

import warnings
from functools import wraps

from . import optimization

__all__ = []


_msg = ("DeprecationWarning: `dask.optimize.{0} has moved to "
        "`dask.optimization.{0}`, please update imports accordingly")


def wrap_deprecated(method):
    __all__.append(method)

    @wraps(getattr(optimization, method))
    def inner(*args, **kwargs):
        warnings.warn(_msg.format(method))
        return getattr(optimization, method)(*args, **kwargs)

    return inner


cull = wrap_deprecated("cull")
fuse_linear = wrap_deprecated("fuse_linear")
inline = wrap_deprecated("inline")
inline_functions = wrap_deprecated("inline_functions")
functions_of = wrap_deprecated("functions_of")
fuse_selections = wrap_deprecated("fuse_selections")
fuse_getitem = wrap_deprecated("fuse_getitem")
fuse = wrap_deprecated("fuse")
key_split = wrap_deprecated("key_split")
