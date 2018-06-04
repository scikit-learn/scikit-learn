from __future__ import absolute_import, division, print_function

from collections import OrderedDict, Iterator
from functools import partial
from hashlib import md5
from operator import getitem
import inspect
import pickle
import os
import threading
import uuid

from toolz import merge, groupby, curry, identity
from toolz.functoolz import Compose

from .compatibility import long, unicode
from .context import _globals, thread_state
from .core import flatten, quote, get as simple_get
from .hashing import hash_buffer_hex
from .utils import Dispatch, ensure_dict


__all__ = ("DaskMethodsMixin",
           "is_dask_collection",
           "compute", "persist", "optimize", "visualize",
           "tokenize", "normalize_token")


def is_dask_collection(x):
    """Returns ``True`` if ``x`` is a dask collection"""
    try:
        return x.__dask_graph__() is not None
    except (AttributeError, TypeError):
        return False


class DaskMethodsMixin(object):
    """A mixin adding standard dask collection methods"""
    __slots__ = ()

    def visualize(self, filename='mydask', format=None, optimize_graph=False,
                  **kwargs):
        """Render the computation of this object's task graph using graphviz.

        Requires ``graphviz`` to be installed.

        Parameters
        ----------
        filename : str or None, optional
            The name (without an extension) of the file to write to disk.  If
            `filename` is None, no file will be written, and we communicate
            with dot using only pipes.
        format : {'png', 'pdf', 'dot', 'svg', 'jpeg', 'jpg'}, optional
            Format in which to write output file.  Default is 'png'.
        optimize_graph : bool, optional
            If True, the graph is optimized before rendering.  Otherwise,
            the graph is displayed as is. Default is False.
        color: {None, 'order'}, optional
            Options to color nodes.  Provide ``cmap=`` keyword for additional
            colormap
        **kwargs
           Additional keyword arguments to forward to ``to_graphviz``.

        Examples
        --------
        >>> x.visualize(filename='dask.pdf')  # doctest: +SKIP
        >>> x.visualize(filename='dask.pdf', color='order')  # doctest: +SKIP

        Returns
        -------
        result : IPython.diplay.Image, IPython.display.SVG, or None
            See dask.dot.dot_graph for more information.

        See Also
        --------
        dask.base.visualize
        dask.dot.dot_graph

        Notes
        -----
        For more information on optimization see here:

        http://dask.pydata.org/en/latest/optimize.html
        """
        return visualize(self, filename=filename, format=format,
                         optimize_graph=optimize_graph, **kwargs)

    def persist(self, **kwargs):
        """Persist this dask collection into memory

        This turns a lazy Dask collection into a Dask collection with the same
        metadata, but now with the results fully computed or actively computing
        in the background.

        The action of function differs significantly depending on the active
        task scheduler.  If the task scheduler supports asynchronous computing,
        such as is the case of the dask.distributed scheduler, then persist
        will return *immediately* and the return value's task graph will
        contain Dask Future objects.  However if the task scheduler only
        supports blocking computation then the call to persist will *block*
        and the return value's task graph will contain concrete Python results.

        This function is particularly useful when using distributed systems,
        because the results will be kept in distributed memory, rather than
        returned to the local process as with compute.

        Parameters
        ----------
        get : callable, optional
            A scheduler ``get`` function to use. If not provided, the default
            is to check the global settings first, and then fall back to
            the collection defaults.
        optimize_graph : bool, optional
            If True [default], the graph is optimized before computation.
            Otherwise the graph is run as is. This can be useful for debugging.
        **kwargs
            Extra keywords to forward to the scheduler ``get`` function.

        Returns
        -------
        New dask collections backed by in-memory data

        See Also
        --------
        dask.base.persist
        """
        (result,) = persist(self, traverse=False, **kwargs)
        return result

    def compute(self, **kwargs):
        """Compute this dask collection

        This turns a lazy Dask collection into its in-memory equivalent.
        For example a Dask.array turns into a  :func:`numpy.array` and a Dask.dataframe
        turns into a Pandas dataframe.  The entire dataset must fit into memory
        before calling this operation.

        Parameters
        ----------
        get : callable, optional
            A scheduler ``get`` function to use. If not provided, the default
            is to check the global settings first, and then fall back to
            the collection defaults.
        optimize_graph : bool, optional
            If True [default], the graph is optimized before computation.
            Otherwise the graph is run as is. This can be useful for debugging.
        kwargs
            Extra keywords to forward to the scheduler ``get`` function.

        See Also
        --------
        dask.base.compute
        """
        (result,) = compute(self, traverse=False, **kwargs)
        return result


def compute_as_if_collection(cls, dsk, keys, get=None, **kwargs):
    """Compute a graph as if it were of type cls.

    Allows for applying the same optimizations and default scheduler."""
    get = get or _globals['get'] or cls.__dask_scheduler__
    dsk2 = optimization_function(cls)(ensure_dict(dsk), keys, **kwargs)
    return get(dsk2, keys, **kwargs)


def dont_optimize(dsk, keys, **kwargs):
    return dsk


def optimization_function(x):
    return getattr(x, '__dask_optimize__', dont_optimize)


def collections_to_dsk(collections, optimize_graph=True, **kwargs):
    """
    Convert many collections into a single dask graph, after optimization
    """
    optimizations = (kwargs.pop('optimizations', None) or
                     _globals.get('optimizations', []))

    if optimize_graph:
        groups = groupby(optimization_function, collections)
        groups = {opt: _extract_graph_and_keys(val)
                  for opt, val in groups.items()}

        for opt in optimizations:
            groups = {k: (opt(dsk, keys), keys)
                      for k, (dsk, keys) in groups.items()}

        dsk = merge(*(opt(dsk, keys, **kwargs)
                      for opt, (dsk, keys) in groups.items()))
    else:
        dsk, _ = _extract_graph_and_keys(collections)

    return dsk


def _extract_graph_and_keys(vals):
    """Given a list of dask vals, return a single graph and a list of keys such
    that ``get(dsk, keys)`` is equivalent to ``[v.compute() v in vals]``."""
    dsk = {}
    keys = []
    for v in vals:
        d = v.__dask_graph__()
        if hasattr(d, 'dicts'):
            for dd in d.dicts.values():
                dsk.update(dd)
        else:
            dsk.update(d)
        keys.append(v.__dask_keys__())

    return dsk, keys


def unpack_collections(*args, **kwargs):
    """Extract collections in preparation for compute/persist/etc...

    Intended use is to find all collections in a set of (possibly nested)
    python objects, do something to them (compute, etc...), then repackage them
    in equivalent python objects.

    Parameters
    ----------
    *args
        Any number of objects. If it is a dask collection, it's extracted and
        added to the list of collections returned. By default, python builtin
        collections are also traversed to look for dask collections (for more
        information see the ``traverse`` keyword).
    traverse : bool, optional
        If True (default), builtin python collections are traversed looking for
        any dask collections they might contain.

    Returns
    -------
    collections : list
        A list of all dask collections contained in ``args``
    repack : callable
        A function to call on the transformed collections to repackage them as
        they were in the original ``args``.
    """
    traverse = kwargs.pop('traverse', True)

    collections = []
    repack_dsk = {}

    def _unpack(expr):
        if is_dask_collection(expr):
            tok = tokenize(expr)
            if tok not in repack_dsk:
                repack_dsk[tok] = (getitem, 'collections', len(collections))
                collections.append(expr)
            return tok

        tok = uuid.uuid4().hex
        if not traverse:
            tsk = quote(expr)
        else:
            # Treat iterators like lists
            typ = list if isinstance(expr, Iterator) else type(expr)

            if typ in (list, tuple, set):
                tsk = (typ, [_unpack(i) for i in expr])
            elif typ is dict:
                tsk = (dict, [[_unpack(k), _unpack(v)]
                              for k, v in expr.items()])
            else:
                return expr

        repack_dsk[tok] = tsk
        return tok

    out = uuid.uuid4().hex
    repack_dsk[out] = (tuple, [_unpack(i) for i in args])

    def repack(results):
        dsk = repack_dsk.copy()
        dsk['collections'] = quote(results)
        return simple_get(dsk, out)

    return collections, repack


def optimize(*args, **kwargs):
    """Optimize several dask collections at once.

    Returns equivalent dask collections that all share the same merged and
    optimized underlying graph. This can be useful if converting multiple
    collections to delayed objects, or to manually apply the optimizations at
    strategic points.

    Note that in most cases you shouldn't need to call this method directly.

    Parameters
    ----------
    *args : objects
        Any number of objects. If a dask object, its graph is optimized and
        merged with all those of all other dask objects before returning an
        equivalent dask collection. Non-dask arguments are passed through
        unchanged.
    traverse : bool, optional
        By default dask traverses builtin python collections looking for dask
        objects passed to ``optimize``. For large collections this can be
        expensive. If none of the arguments contain any dask objects, set
        ``traverse=False`` to avoid doing this traversal.
    optimizations : list of callables, optional
        Additional optimization passes to perform.
    **kwargs
        Extra keyword arguments to forward to the optimization passes.

    Examples
    --------
    >>> import dask.array as da
    >>> a = da.arange(10, chunks=2).sum()
    >>> b = da.arange(10, chunks=2).mean()
    >>> a2, b2 = optimize(a, b)

    >>> a2.compute() == a.compute()
    True
    >>> b2.compute() == b.compute()
    True
    """
    collections, repack = unpack_collections(*args, **kwargs)
    if not collections:
        return args

    dsk = collections_to_dsk(collections, **kwargs)
    postpersists = [a.__dask_postpersist__() if is_dask_collection(a)
                    else (None, a) for a in args]

    keys, postpersists = [], []
    for a in collections:
        keys.extend(flatten(a.__dask_keys__()))
        postpersists.append(a.__dask_postpersist__())

    return repack([r(dsk, *s) for r, s in postpersists])


# TODO: remove after deprecation cycle of `dask.optimize` module completes
from . import optimize as _deprecated_optimize
for _m in _deprecated_optimize.__all__:
    setattr(optimize, _m, getattr(_deprecated_optimize, _m))


def compute(*args, **kwargs):
    """Compute several dask collections at once.

    Parameters
    ----------
    args : object
        Any number of objects. If it is a dask object, it's computed and the
        result is returned. By default, python builtin collections are also
        traversed to look for dask objects (for more information see the
        ``traverse`` keyword). Non-dask arguments are passed through unchanged.
    traverse : bool, optional
        By default dask traverses builtin python collections looking for dask
        objects passed to ``compute``. For large collections this can be
        expensive. If none of the arguments contain any dask objects, set
        ``traverse=False`` to avoid doing this traversal.
    get : callable, optional
        A scheduler ``get`` function to use. If not provided, the default is
        to check the global settings first, and then fall back to defaults for
        the collections.
    optimize_graph : bool, optional
        If True [default], the optimizations for each collection are applied
        before computation. Otherwise the graph is run as is. This can be
        useful for debugging.
    kwargs
        Extra keywords to forward to the scheduler ``get`` function.

    Examples
    --------
    >>> import dask.array as da
    >>> a = da.arange(10, chunks=2).sum()
    >>> b = da.arange(10, chunks=2).mean()
    >>> compute(a, b)
    (45, 4.5)

    By default, dask objects inside python collections will also be computed:

    >>> compute({'a': a, 'b': b, 'c': 1})  # doctest: +SKIP
    ({'a': 45, 'b': 4.5, 'c': 1},)
    """
    traverse = kwargs.pop('traverse', True)
    optimize_graph = kwargs.pop('optimize_graph', True)
    get = kwargs.pop('get', None) or _globals['get']

    collections, repack = unpack_collections(*args, traverse=traverse)
    if not collections:
        return args

    if get is None and getattr(thread_state, 'key', False):
        from distributed.worker import get_worker
        get = get_worker().client.get

    if not get:
        get = collections[0].__dask_scheduler__
        if not all(a.__dask_scheduler__ == get for a in collections):
            raise ValueError("Compute called on multiple collections with "
                             "differing default schedulers. Please specify a "
                             "scheduler `get` function using either "
                             "the `get` kwarg or globally with `set_options`.")

    dsk = collections_to_dsk(collections, optimize_graph, **kwargs)
    keys = [x.__dask_keys__() for x in collections]
    postcomputes = [x.__dask_postcompute__() for x in collections]
    results = get(dsk, keys, **kwargs)
    return repack([f(r, *a) for r, (f, a) in zip(results, postcomputes)])


def visualize(*args, **kwargs):
    """
    Visualize several dask graphs at once.

    Requires ``graphviz`` to be installed. All options that are not the dask
    graph(s) should be passed as keyword arguments.

    Parameters
    ----------
    dsk : dict(s) or collection(s)
        The dask graph(s) to visualize.
    filename : str or None, optional
        The name (without an extension) of the file to write to disk.  If
        `filename` is None, no file will be written, and we communicate
        with dot using only pipes.
    format : {'png', 'pdf', 'dot', 'svg', 'jpeg', 'jpg'}, optional
        Format in which to write output file.  Default is 'png'.
    optimize_graph : bool, optional
        If True, the graph is optimized before rendering.  Otherwise,
        the graph is displayed as is. Default is False.
    color: {None, 'order'}, optional
        Options to color nodes.  Provide ``cmap=`` keyword for additional
        colormap
    **kwargs
       Additional keyword arguments to forward to ``to_graphviz``.

    Examples
    --------
    >>> x.visualize(filename='dask.pdf')  # doctest: +SKIP
    >>> x.visualize(filename='dask.pdf', color='order')  # doctest: +SKIP

    Returns
    -------
    result : IPython.diplay.Image, IPython.display.SVG, or None
        See dask.dot.dot_graph for more information.

    See Also
    --------
    dask.dot.dot_graph

    Notes
    -----
    For more information on optimization see here:

    http://dask.pydata.org/en/latest/optimize.html
    """
    from dask.dot import dot_graph

    filename = kwargs.pop('filename', 'mydask')
    optimize_graph = kwargs.pop('optimize_graph', False)

    dsks = [arg for arg in args if isinstance(arg, dict)]
    args = [arg for arg in args if is_dask_collection(arg)]

    dsk = collections_to_dsk(args, optimize_graph=optimize_graph)
    for d in dsks:
        dsk.update(d)

    color = kwargs.get('color')

    if color == 'order':
        from .order import order
        import matplotlib.pyplot as plt
        o = order(dsk)
        try:
            cmap = kwargs.pop('cmap')
        except KeyError:
            cmap = plt.cm.RdBu
        if isinstance(cmap, str):
            import matplotlib.pyplot as plt
            cmap = getattr(plt.cm, cmap)
        mx = max(o.values()) + 1
        colors = {k: _colorize(cmap(v / mx, bytes=True)) for k, v in o.items()}

        kwargs['function_attributes'] = {k: {'color': v, 'label': str(o[k])}
                                         for k, v in colors.items()}
        kwargs['data_attributes'] = {k: {'color': v} for k, v in colors.items()}
    elif color:
        raise NotImplementedError("Unknown value color=%s" % color)

    return dot_graph(dsk, filename=filename, **kwargs)


def persist(*args, **kwargs):
    """ Persist multiple Dask collections into memory

    This turns lazy Dask collections into Dask collections with the same
    metadata, but now with their results fully computed or actively computing
    in the background.

    For example a lazy dask.array built up from many lazy calls will now be a
    dask.array of the same shape, dtype, chunks, etc., but now with all of
    those previously lazy tasks either computed in memory as many small :class:`numpy.array`
    (in the single-machine case) or asynchronously running in the
    background on a cluster (in the distributed case).

    This function operates differently if a ``dask.distributed.Client`` exists
    and is connected to a distributed scheduler.  In this case this function
    will return as soon as the task graph has been submitted to the cluster,
    but before the computations have completed.  Computations will continue
    asynchronously in the background.  When using this function with the single
    machine scheduler it blocks until the computations have finished.

    When using Dask on a single machine you should ensure that the dataset fits
    entirely within memory.

    Examples
    --------
    >>> df = dd.read_csv('/path/to/*.csv')  # doctest: +SKIP
    >>> df = df[df.name == 'Alice']  # doctest: +SKIP
    >>> df['in-debt'] = df.balance < 0  # doctest: +SKIP
    >>> df = df.persist()  # triggers computation  # doctest: +SKIP

    >>> df.value().min()  # future computations are now fast  # doctest: +SKIP
    -10
    >>> df.value().max()  # doctest: +SKIP
    100

    >>> from dask import persist  # use persist function on multiple collections
    >>> a, b = persist(a, b)  # doctest: +SKIP

    Parameters
    ----------
    *args: Dask collections
    get : callable, optional
        A scheduler ``get`` function to use. If not provided, the default
        is to check the global settings first, and then fall back to
        the collection defaults.
    traverse : bool, optional
        By default dask traverses builtin python collections looking for dask
        objects passed to ``persist``. For large collections this can be
        expensive. If none of the arguments contain any dask objects, set
        ``traverse=False`` to avoid doing this traversal.
    optimize_graph : bool, optional
        If True [default], the graph is optimized before computation.
        Otherwise the graph is run as is. This can be useful for debugging.
    **kwargs
        Extra keywords to forward to the scheduler ``get`` function.

    Returns
    -------
    New dask collections backed by in-memory data
    """
    traverse = kwargs.pop('traverse', True)
    optimize_graph = kwargs.pop('optimize_graph', True)
    get = kwargs.pop('get', None) or _globals['get']

    collections, repack = unpack_collections(*args, traverse=traverse)
    if not collections:
        return args

    if get is None and getattr(thread_state, 'key', False):
        from distributed.worker import get_worker
        get = get_worker().client.get

    if inspect.ismethod(get):
        try:
            from distributed.client import default_client
        except ImportError:
            pass
        else:
            try:
                client = default_client()
            except ValueError:
                pass
            else:
                if client.get == get:
                    results = client.persist(collections,
                                             optimize_graph=optimize_graph,
                                             **kwargs)
                    return repack(results)

    if not get:
        get = collections[0].__dask_scheduler__
        if not all(a.__dask_scheduler__ == get for a in collections):
            raise ValueError("Persist called on multiple collections with "
                             "differing default schedulers. Please specify a "
                             "scheduler `get` function using either "
                             "the `get` kwarg or globally with `set_options`.")

    dsk = collections_to_dsk(collections, optimize_graph, **kwargs)
    keys, postpersists = [], []
    for a in collections:
        a_keys = list(flatten(a.__dask_keys__()))
        rebuild, state = a.__dask_postpersist__()
        keys.extend(a_keys)
        postpersists.append((rebuild, a_keys, state))

    results = get(dsk, keys, **kwargs)
    d = dict(zip(keys, results))
    results2 = [r({k: d[k] for k in ks}, *s) for r, ks, s in postpersists]
    return repack(results2)


############
# Tokenize #
############

def tokenize(*args, **kwargs):
    """ Deterministic token

    >>> tokenize([1, 2, '3'])
    '7d6a880cd9ec03506eee6973ff551339'

    >>> tokenize('Hello') == tokenize('Hello')
    True
    """
    if kwargs:
        args = args + (kwargs,)
    return md5(str(tuple(map(normalize_token, args))).encode()).hexdigest()


normalize_token = Dispatch()
normalize_token.register((int, long, float, str, unicode, bytes, type(None),
                          type, slice, complex, type(Ellipsis)),
                         identity)


@normalize_token.register(dict)
def normalize_dict(d):
    return normalize_token(sorted(d.items(), key=str))


@normalize_token.register(OrderedDict)
def normalize_ordered_dict(d):
    return type(d).__name__, normalize_token(list(d.items()))


@normalize_token.register(set)
def normalize_set(s):
    return normalize_token(sorted(s, key=str))


@normalize_token.register((tuple, list))
def normalize_seq(seq):
    return type(seq).__name__, list(map(normalize_token, seq))


@normalize_token.register(object)
def normalize_object(o):
    method = getattr(o, '__dask_tokenize__', None)
    if method is not None:
        return method()
    return normalize_function(o) if callable(o) else uuid.uuid4().hex


function_cache = {}
function_cache_lock = threading.Lock()


def normalize_function(func):
    try:
        return function_cache[func]
    except KeyError:
        result = _normalize_function(func)
        if len(function_cache) >= 500:  # clear half of cache if full
            with function_cache_lock:
                if len(function_cache) >= 500:
                    for k in list(function_cache)[::2]:
                        del function_cache[k]
        function_cache[func] = result
        return result
    except TypeError:  # not hashable
        return _normalize_function(func)


def _normalize_function(func):
    if isinstance(func, curry):
        func = func._partial
    if isinstance(func, Compose):
        first = getattr(func, 'first', None)
        funcs = reversed((first,) + func.funcs) if first else func.funcs
        return tuple(normalize_function(f) for f in funcs)
    elif isinstance(func, partial):
        args = tuple(normalize_token(i) for i in func.args)
        if func.keywords:
            kws = tuple((k, normalize_token(v))
                        for k, v in sorted(func.keywords.items()))
        else:
            kws = None
        return (normalize_function(func.func), args, kws)
    else:
        try:
            result = pickle.dumps(func, protocol=0)
            if b'__main__' not in result:  # abort on dynamic functions
                return result
        except Exception:
            pass
        try:
            import cloudpickle
            return cloudpickle.dumps(func, protocol=0)
        except Exception:
            return str(func)


@normalize_token.register_lazy("pandas")
def register_pandas():
    import pandas as pd

    @normalize_token.register(pd.Index)
    def normalize_index(ind):
        return [ind.name, normalize_token(ind.values)]

    @normalize_token.register(pd.Categorical)
    def normalize_categorical(cat):
        return [normalize_token(cat.codes),
                normalize_token(cat.categories),
                cat.ordered]

    @normalize_token.register(pd.Series)
    def normalize_series(s):
        return [s.name, s.dtype,
                normalize_token(s._data.blocks[0].values),
                normalize_token(s.index)]

    @normalize_token.register(pd.DataFrame)
    def normalize_dataframe(df):
        data = [block.values for block in df._data.blocks]
        data += [df.columns, df.index]
        return list(map(normalize_token, data))


@normalize_token.register_lazy("numpy")
def register_numpy():
    import numpy as np

    @normalize_token.register(np.ndarray)
    def normalize_array(x):
        if not x.shape:
            return (str(x), x.dtype)
        if hasattr(x, 'mode') and getattr(x, 'filename', None):
            if hasattr(x.base, 'ctypes'):
                offset = (x.ctypes.get_as_parameter().value -
                          x.base.ctypes.get_as_parameter().value)
            else:
                offset = 0  # root memmap's have mmap object as base
            return (x.filename, os.path.getmtime(x.filename), x.dtype,
                    x.shape, x.strides, offset)
        if x.dtype.hasobject:
            try:
                data = hash_buffer_hex('-'.join(x.flat).encode('utf-8'))
            except TypeError:
                data = hash_buffer_hex(b'-'.join([unicode(item).encode('utf-8') for item in
                                                  x.flat]))
        else:
            try:
                data = hash_buffer_hex(x.ravel(order='K').view('i1'))
            except (BufferError, AttributeError, ValueError):
                data = hash_buffer_hex(x.copy().ravel(order='K').view('i1'))
        return (data, x.dtype, x.shape, x.strides)

    @normalize_token.register(np.matrix)
    def normalize_matrix(x):
        return type(x).__name__, normalize_array(x.view(type=np.ndarray))

    normalize_token.register(np.dtype, repr)
    normalize_token.register(np.generic, repr)

    @normalize_token.register(np.ufunc)
    def normalize_ufunc(x):
        try:
            name = x.__name__
            if getattr(np, name) is x:
                return 'np.' + name
        except AttributeError:
            return normalize_function(x)


@normalize_token.register_lazy("scipy")
def register_scipy():
    import scipy.sparse as sp

    def normalize_sparse_matrix(x, attrs):
        return type(x).__name__, normalize_seq((normalize_token(getattr(x, key))
                                                for key in attrs))

    for cls, attrs in [(sp.dia_matrix, ('data', 'offsets', 'shape')),
                       (sp.bsr_matrix, ('data', 'indices', 'indptr',
                                        'blocksize', 'shape')),
                       (sp.coo_matrix, ('data', 'row', 'col', 'shape')),
                       (sp.csr_matrix, ('data', 'indices', 'indptr', 'shape')),
                       (sp.csc_matrix, ('data', 'indices', 'indptr', 'shape')),
                       (sp.lil_matrix, ('data', 'rows', 'shape'))]:
        normalize_token.register(cls,
                                 partial(normalize_sparse_matrix, attrs=attrs))

    @normalize_token.register(sp.dok_matrix)
    def normalize_dok_matrix(x):
        return type(x).__name__, normalize_token(sorted(x.items()))


def _colorize(t):
    """ Convert (r, g, b) triple to "#RRGGBB" string

    For use with ``visualize(color=...)``

    Examples
    --------
    >>> _colorize((255, 255, 255))
    '#FFFFFF'
    >>> _colorize((0, 32, 128))
    '#002080'
    """
    t = t[:3]
    i = sum(v * 256 ** (len(t) - i - 1) for i, v in enumerate(t))
    h = hex(int(i))[2:].upper()
    h = '0' * (6 - len(h)) + h
    return "#" + h
