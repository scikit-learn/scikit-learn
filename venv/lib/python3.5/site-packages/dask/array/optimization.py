from __future__ import absolute_import, division, print_function

from operator import getitem

import numpy as np

from .core import getter, getter_nofancy, getter_inline
from ..compatibility import zip_longest
from ..core import flatten, reverse_dict
from ..optimization import cull, fuse, inline_functions
from ..utils import ensure_dict


# All get* functions the optimizations know about
GETTERS = (getter, getter_nofancy, getter_inline, getitem)
# These get* functions aren't ever completely removed from the graph,
# even if the index should be a no-op by numpy semantics. Some array-like's
# don't completely follow semantics, making indexing always necessary.
GETNOREMOVE = (getter, getter_nofancy)


def optimize(dsk, keys, fuse_keys=None, fast_functions=None,
             inline_functions_fast_functions=(getter_inline,), rename_fused_keys=True,
             **kwargs):
    """ Optimize dask for array computation

    1.  Cull tasks not necessary to evaluate keys
    2.  Remove full slicing, e.g. x[:]
    3.  Inline fast functions like getitem and np.transpose
    """
    dsk = ensure_dict(dsk)
    keys = list(flatten(keys))
    if fast_functions is not None:
        inline_functions_fast_functions = fast_functions

    dsk2, dependencies = cull(dsk, keys)
    hold = hold_keys(dsk2, dependencies)

    dsk3, dependencies = fuse(dsk2, hold + keys + (fuse_keys or []),
                              dependencies, rename_keys=rename_fused_keys)
    if inline_functions_fast_functions:
        dsk4 = inline_functions(dsk3, keys, dependencies=dependencies,
                                fast_functions=inline_functions_fast_functions)
    else:
        dsk4 = dsk3
    dsk5 = optimize_slices(dsk4)

    return dsk5


def hold_keys(dsk, dependencies):
    """ Find keys to avoid fusion

    We don't want to fuse data present in the graph because it is easier to
    serialize as a raw value.

    We don't want to fuse chains after getitem/GETTERS because we want to
    move around only small pieces of data, rather than the underlying arrays.
    """
    dependents = reverse_dict(dependencies)
    data = {k for k, v in dsk.items() if type(v) not in (tuple, str)}

    hold_keys = list(data)
    for dat in data:
        deps = dependents[dat]
        for dep in deps:
            task = dsk[dep]
            # If the task is a get* function, we walk up the chain, and stop
            # when there's either more than one dependent, or the dependent is
            # no longer a get* function or an alias. We then add the final
            # key to the list of keys not to fuse.
            if type(task) is tuple and task and task[0] in GETTERS:
                try:
                    while len(dependents[dep]) == 1:
                        new_dep = next(iter(dependents[dep]))
                        new_task = dsk[new_dep]
                        # If the task is a get* or an alias, continue up the
                        # linear chain
                        if new_task[0] in GETTERS or new_task in dsk:
                            dep = new_dep
                        else:
                            break
                except (IndexError, TypeError):
                    pass
                hold_keys.append(dep)
    return hold_keys


def optimize_slices(dsk):
    """ Optimize slices

    1.  Fuse repeated slices, like x[5:][2:6] -> x[7:11]
    2.  Remove full slices, like         x[:] -> x

    See also:
        fuse_slice_dict
    """
    fancy_ind_types = (list, np.ndarray)
    dsk = dsk.copy()
    for k, v in dsk.items():
        if type(v) is tuple and v[0] in GETTERS and len(v) in (3, 5):
            if len(v) == 3:
                get, a, a_index = v
                # getter defaults to asarray=True, getitem is semantically False
                a_asarray = get is not getitem
                a_lock = None
            else:
                get, a, a_index, a_asarray, a_lock = v
            while type(a) is tuple and a[0] in GETTERS and len(a) in (3, 5):
                if len(a) == 3:
                    f2, b, b_index = a
                    b_asarray = f2 is not getitem
                    b_lock = None
                else:
                    f2, b, b_index, b_asarray, b_lock = a

                if a_lock and a_lock is not b_lock:
                    break
                if (type(a_index) is tuple) != (type(b_index) is tuple):
                    break
                if type(a_index) is tuple:
                    indices = b_index + a_index
                    if (len(a_index) != len(b_index) and
                            any(i is None for i in indices)):
                        break
                    if (f2 is getter_nofancy and
                            any(isinstance(i, fancy_ind_types) for i in indices)):
                        break
                elif (f2 is getter_nofancy and
                        (type(a_index) in fancy_ind_types or
                         type(b_index) in fancy_ind_types)):
                    break
                try:
                    c_index = fuse_slice(b_index, a_index)
                    # rely on fact that nested gets never decrease in
                    # strictness e.g. `(getter_nofancy, (getter, ...))` never
                    # happens
                    get = getter if f2 is getter_inline else f2
                except NotImplementedError:
                    break
                a, a_index, a_lock = b, c_index, b_lock
                a_asarray |= b_asarray

            # Skip the get call if not from from_array and nothing to do
            if (get not in GETNOREMOVE and
                ((type(a_index) is slice and not a_index.start and
                  a_index.stop is None and a_index.step is None) or
                 (type(a_index) is tuple and
                  all(type(s) is slice and not s.start and s.stop is None and
                      s.step is None for s in a_index)))):
                dsk[k] = a
            elif get is getitem or (a_asarray and not a_lock):
                # default settings are fine, drop the extra parameters Since we
                # always fallback to inner `get` functions, `get is getitem`
                # can only occur if all gets are getitem, meaning all
                # parameters must be getitem defaults.
                dsk[k] = (get, a, a_index)
            else:
                dsk[k] = (get, a, a_index, a_asarray, a_lock)

    return dsk


def normalize_slice(s):
    """ Replace Nones in slices with integers

    >>> normalize_slice(slice(None, None, None))
    slice(0, None, 1)
    """
    start, stop, step = s.start, s.stop, s.step
    if start is None:
        start = 0
    if step is None:
        step = 1
    if start < 0 or step < 0 or stop is not None and stop < 0:
        raise NotImplementedError()
    return slice(start, stop, step)


def check_for_nonfusible_fancy_indexing(fancy, normal):
    # Check for fancy indexing and normal indexing, where the fancy
    # indexed dimensions != normal indexed dimensions with integers. E.g.:
    # disallow things like:
    # x[:, [1, 2], :][0, :, :] -> x[0, [1, 2], :] or
    # x[0, :, :][:, [1, 2], :] -> x[0, [1, 2], :]
    for f, n in zip_longest(fancy, normal, fillvalue=slice(None)):
        if type(f) is not list and isinstance(n, int):
            raise NotImplementedError("Can't handle normal indexing with "
                                      "integers and fancy indexing if the "
                                      "integers and fancy indices don't "
                                      "align with the same dimensions.")


def fuse_slice(a, b):
    """ Fuse stacked slices together

    Fuse a pair of repeated slices into a single slice:

    >>> fuse_slice(slice(1000, 2000), slice(10, 15))
    slice(1010, 1015, None)

    This also works for tuples of slices

    >>> fuse_slice((slice(100, 200), slice(100, 200, 10)),
    ...            (slice(10, 15), [5, 2]))
    (slice(110, 115, None), [150, 120])

    And a variety of other interesting cases

    >>> fuse_slice(slice(1000, 2000), 10)  # integers
    1010

    >>> fuse_slice(slice(1000, 2000, 5), slice(10, 20, 2))
    slice(1050, 1100, 10)

    >>> fuse_slice(slice(1000, 2000, 5), [1, 2, 3])  # lists
    [1005, 1010, 1015]

    >>> fuse_slice(None, slice(None, None))  # doctest: +SKIP
    None
    """
    # None only works if the second side is a full slice
    if a is None and isinstance(b, slice) and b == slice(None, None):
        return None

    # Replace None with 0 and one in start and step
    if isinstance(a, slice):
        a = normalize_slice(a)
    if isinstance(b, slice):
        b = normalize_slice(b)

    if isinstance(a, slice) and isinstance(b, int):
        if b < 0:
            raise NotImplementedError()
        return a.start + b * a.step

    if isinstance(a, slice) and isinstance(b, slice):
        start = a.start + a.step * b.start
        if b.stop is not None:
            stop = a.start + a.step * b.stop
        else:
            stop = None
        if a.stop is not None:
            if stop is not None:
                stop = min(a.stop, stop)
            else:
                stop = a.stop
        step = a.step * b.step
        if step == 1:
            step = None
        return slice(start, stop, step)

    if isinstance(b, list):
        return [fuse_slice(a, bb) for bb in b]
    if isinstance(a, list) and isinstance(b, (int, slice)):
        return a[b]

    if isinstance(a, tuple) and not isinstance(b, tuple):
        b = (b,)

    # If given two tuples walk through both, being mindful of uneven sizes
    # and newaxes
    if isinstance(a, tuple) and isinstance(b, tuple):

        # Check for non-fusible cases with fancy-indexing
        a_has_lists = any(isinstance(item, list) for item in a)
        b_has_lists = any(isinstance(item, list) for item in b)
        if a_has_lists and b_has_lists:
            raise NotImplementedError("Can't handle multiple list indexing")
        elif a_has_lists:
            check_for_nonfusible_fancy_indexing(a, b)
        elif b_has_lists:
            check_for_nonfusible_fancy_indexing(b, a)

        j = 0
        result = list()
        for i in range(len(a)):
            #  axis ceased to exist  or we're out of b
            if isinstance(a[i], int) or j == len(b):
                result.append(a[i])
                continue
            while b[j] is None:  # insert any Nones on the rhs
                result.append(None)
                j += 1
            result.append(fuse_slice(a[i], b[j]))  # Common case
            j += 1
        while j < len(b):  # anything leftover on the right?
            result.append(b[j])
            j += 1
        return tuple(result)
    raise NotImplementedError()
