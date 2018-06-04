from __future__ import absolute_import, division, print_function

from functools import partial
from itertools import product

import numpy as np

try:
    from cytoolz import curry
except ImportError:
    from toolz import curry

from ..base import tokenize
from .core import Array, normalize_chunks


def wrap_func_shape_as_first_arg(func, *args, **kwargs):
    """
    Transform np creation function into blocked version
    """
    if 'shape' not in kwargs:
        shape, args = args[0], args[1:]
    else:
        shape = kwargs.pop('shape')

    if not isinstance(shape, (tuple, list)):
        shape = (shape,)

    chunks = kwargs.pop('chunks', None)
    chunks = normalize_chunks(chunks, shape)
    name = kwargs.pop('name', None)

    dtype = kwargs.pop('dtype', None)
    if dtype is None:
        dtype = func(shape, *args, **kwargs).dtype

    name = name or 'wrapped-' + tokenize(func, shape, chunks, dtype, args, kwargs)

    keys = product([name], *[range(len(bd)) for bd in chunks])
    shapes = product(*chunks)
    func = partial(func, dtype=dtype, **kwargs)
    vals = ((func,) + (s,) + args for s in shapes)

    dsk = dict(zip(keys, vals))
    return Array(dsk, name, chunks, dtype=dtype)


@curry
def wrap(wrap_func, func, **kwargs):
    f = partial(wrap_func, func, **kwargs)
    template = """
    Blocked variant of %(name)s

    Follows the signature of %(name)s exactly except that it also requires a
    keyword argument chunks=(...)

    Original signature follows below.
    """
    if func.__doc__ is not None:
        f.__doc__ = template % {'name': func.__name__} + func.__doc__
        f.__name__ = 'blocked_' + func.__name__
    return f


w = wrap(wrap_func_shape_as_first_arg)

ones = w(np.ones, dtype='f8')
zeros = w(np.zeros, dtype='f8')
empty = w(np.empty, dtype='f8')
full = w(np.full)
