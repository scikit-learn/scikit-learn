from __future__ import absolute_import, division, print_function

import io
import itertools
import math
import types
import uuid
import warnings
from collections import Iterable, Iterator, defaultdict
from distutils.version import LooseVersion
from functools import wraps, partial
from operator import getitem
from random import Random

from toolz import (merge, take, reduce, valmap, map, partition_all, filter,
                   remove, compose, curry, first, second, accumulate, peek)
from toolz.compatibility import iteritems, zip
import toolz
_implement_accumulate = LooseVersion(toolz.__version__) > '0.7.4'
try:
    import cytoolz
    from cytoolz import (frequencies, merge_with, join, reduceby,
                         count, pluck, groupby, topk)
    if LooseVersion(cytoolz.__version__) > '0.7.3':
        from cytoolz import accumulate  # noqa: F811
        _implement_accumulate = True
except ImportError:
    from toolz import (frequencies, merge_with, join, reduceby,
                       count, pluck, groupby, topk)

from ..base import tokenize, dont_optimize, is_dask_collection, DaskMethodsMixin
from ..bytes import open_files
from ..compatibility import apply, urlopen
from ..context import _globals, globalmethod
from ..core import quote, istask, get_dependencies, reverse_dict
from ..delayed import Delayed
from ..multiprocessing import get as mpget
from ..optimization import fuse, cull, inline
from ..utils import (system_encoding, takes_multiple_arguments, funcname,
                     digit, insert, ensure_dict, ensure_bytes, ensure_unicode)


no_default = '__no__default__'
no_result = type('no_result', (object,),
                 {'__slots__': (),
                  '__reduce__': lambda self: 'no_result'})


def lazify_task(task, start=True):
    """
    Given a task, remove unnecessary calls to ``list`` and ``reify``.

    This traverses tasks and small lists.  We choose not to traverse down lists
    of size >= 50 because it is unlikely that sequences this long contain other
    sequences in practice.

    Examples
    --------
    >>> task = (sum, (list, (map, inc, [1, 2, 3])))  # doctest: +SKIP
    >>> lazify_task(task)  # doctest: +SKIP
    (sum, (map, inc, [1, 2, 3]))
    """
    if type(task) is list and len(task) < 50:
        return [lazify_task(arg, False) for arg in task]
    if not istask(task):
        return task
    head, tail = task[0], task[1:]
    if not start and head in (list, reify):
        task = task[1]
        return lazify_task(*tail, start=False)
    else:
        return (head,) + tuple([lazify_task(arg, False) for arg in tail])


def lazify(dsk):
    """
    Remove unnecessary calls to ``list`` in tasks.

    See Also
    --------
    ``dask.bag.core.lazify_task``
    """
    return valmap(lazify_task, dsk)


def inline_singleton_lists(dsk, dependencies=None):
    """ Inline lists that are only used once.

    >>> d = {'b': (list, 'a'),
    ...      'c': (f, 'b', 1)}     # doctest: +SKIP
    >>> inline_singleton_lists(d)  # doctest: +SKIP
    {'c': (f, (list, 'a'), 1)}

    Pairs nicely with lazify afterwards.
    """
    if dependencies is None:
        dependencies = {k: get_dependencies(dsk, task=v)
                        for k, v in dsk.items()}
    dependents = reverse_dict(dependencies)

    keys = [k for k, v in dsk.items()
            if istask(v) and v and v[0] is list and len(dependents[k]) == 1]
    dsk = inline(dsk, keys, inline_constants=False)
    for k in keys:
        del dsk[k]
    return dsk


def optimize(dsk, keys, fuse_keys=None, rename_fused_keys=True, **kwargs):
    """ Optimize a dask from a dask Bag. """
    dsk2, dependencies = cull(dsk, keys)
    dsk3, dependencies = fuse(dsk2, keys + (fuse_keys or []), dependencies,
                              rename_keys=rename_fused_keys)
    dsk4 = inline_singleton_lists(dsk3, dependencies)
    dsk5 = lazify(dsk4)
    return dsk5


def _to_textfiles_chunk(data, lazy_file):
    with lazy_file as f:
        if isinstance(f, io.TextIOWrapper):
            endline = u'\n'
            ensure = ensure_unicode
        else:
            endline = b'\n'
            ensure = ensure_bytes
        started = False
        for d in data:
            if started:
                f.write(endline)
            else:
                started = True
            f.write(ensure(d))


def to_textfiles(b, path, name_function=None, compression='infer',
                 encoding=system_encoding, compute=True, get=None,
                 storage_options=None):
    """ Write dask Bag to disk, one filename per partition, one line per element.

    **Paths**: This will create one file for each partition in your bag. You
    can specify the filenames in a variety of ways.

    Use a globstring

    >>> b.to_textfiles('/path/to/data/*.json.gz')  # doctest: +SKIP

    The * will be replaced by the increasing sequence 1, 2, ...

    ::

        /path/to/data/0.json.gz
        /path/to/data/1.json.gz

    Use a globstring and a ``name_function=`` keyword argument.  The
    name_function function should expect an integer and produce a string.
    Strings produced by name_function must preserve the order of their
    respective partition indices.

    >>> from datetime import date, timedelta
    >>> def name(i):
    ...     return str(date(2015, 1, 1) + i * timedelta(days=1))

    >>> name(0)
    '2015-01-01'
    >>> name(15)
    '2015-01-16'

    >>> b.to_textfiles('/path/to/data/*.json.gz', name_function=name)  # doctest: +SKIP

    ::

        /path/to/data/2015-01-01.json.gz
        /path/to/data/2015-01-02.json.gz
        ...

    You can also provide an explicit list of paths.

    >>> paths = ['/path/to/data/alice.json.gz', '/path/to/data/bob.json.gz', ...]  # doctest: +SKIP
    >>> b.to_textfiles(paths) # doctest: +SKIP

    **Compression**: Filenames with extensions corresponding to known
    compression algorithms (gz, bz2) will be compressed accordingly.

    **Bag Contents**: The bag calling ``to_textfiles`` must be a bag of
    text strings. For example, a bag of dictionaries could be written to
    JSON text files by mapping ``json.dumps`` on to the bag first, and
    then calling ``to_textfiles`` :

    >>> b_dict.map(json.dumps).to_textfiles("/path/to/data/*.json")  # doctest: +SKIP
    """
    mode = 'wb' if encoding is None else 'wt'
    files = open_files(path, compression=compression, mode=mode,
                       encoding=encoding, name_function=name_function,
                       num=b.npartitions, **(storage_options or {}))

    name = 'to-textfiles-' + uuid.uuid4().hex
    dsk = {(name, i): (_to_textfiles_chunk, (b.name, i), f)
           for i, f in enumerate(files)}
    out = type(b)(merge(dsk, b.dask), name, b.npartitions)

    if compute:
        out.compute(get=get)
        return [f.path for f in files]
    else:
        return out.to_delayed()


def finalize(results):
    if not results:
        return results
    if isinstance(results, Iterator):
        results = list(results)
    if isinstance(results[0], Iterable) and not isinstance(results[0], str):
        results = toolz.concat(results)
    if isinstance(results, Iterator):
        results = list(results)
    return results


def finalize_item(results):
    return results[0]


class StringAccessor(object):
    """ String processing functions

    Examples
    --------

    >>> import dask.bag as db
    >>> b = db.from_sequence(['Alice Smith', 'Bob Jones', 'Charlie Smith'])
    >>> list(b.str.lower())
    ['alice smith', 'bob jones', 'charlie smith']

    >>> list(b.str.match('*Smith'))
    ['Alice Smith', 'Charlie Smith']

    >>> list(b.str.split(' '))
    [['Alice', 'Smith'], ['Bob', 'Jones'], ['Charlie', 'Smith']]
    """
    def __init__(self, bag):
        self._bag = bag

    def __dir__(self):
        return sorted(set(dir(type(self)) + dir(str)))

    def _strmap(self, key, *args, **kwargs):
        return self._bag.map(lambda s: getattr(s, key)(*args, **kwargs))

    def __getattr__(self, key):
        try:
            return object.__getattribute__(self, key)
        except AttributeError:
            if key in dir(str):
                func = getattr(str, key)
                return robust_wraps(func)(partial(self._strmap, key))
            else:
                raise

    def match(self, pattern):
        """ Filter strings by those that match a pattern.

        Examples
        --------

        >>> import dask.bag as db
        >>> b = db.from_sequence(['Alice Smith', 'Bob Jones', 'Charlie Smith'])
        >>> list(b.str.match('*Smith'))
        ['Alice Smith', 'Charlie Smith']

        See Also
        --------
        fnmatch.fnmatch
        """
        from fnmatch import fnmatch
        return self._bag.filter(partial(fnmatch, pat=pattern))


def robust_wraps(wrapper):
    """ A weak version of wraps that only copies doc. """
    def _(wrapped):
        wrapped.__doc__ = wrapper.__doc__
        return wrapped
    return _


class Item(DaskMethodsMixin):
    def __init__(self, dsk, key):
        self.dask = dsk
        self.key = key
        self.name = key

    def __dask_graph__(self):
        return self.dask

    def __dask_keys__(self):
        return [self.key]

    def __dask_tokenize__(self):
        return self.key

    __dask_optimize__ = globalmethod(optimize, key='bag_optimize',
                                     falsey=dont_optimize)
    __dask_scheduler__ = staticmethod(mpget)

    def __dask_postcompute__(self):
        return finalize_item, ()

    def __dask_postpersist__(self):
        return Item, (self.key,)

    @staticmethod
    def from_delayed(value):
        """ Create bag item from a dask.delayed value.

        See ``dask.bag.from_delayed`` for details
        """
        from dask.delayed import Delayed, delayed
        if not isinstance(value, Delayed) and hasattr(value, 'key'):
            value = delayed(value)
        assert isinstance(value, Delayed)
        return Item(ensure_dict(value.dask), value.key)

    @property
    def _args(self):
        return (self.dask, self.key)

    def __getstate__(self):
        return self._args

    def __setstate__(self, state):
        self.dask, self.key = state

    def apply(self, func):
        name = 'apply-{0}-{1}'.format(funcname(func), tokenize(self, func))
        dsk = {name: (func, self.key)}
        return Item(merge(self.dask, dsk), name)

    __int__ = __float__ = __complex__ = __bool__ = DaskMethodsMixin.compute

    def to_delayed(self, optimize_graph=True):
        """Convert into a ``dask.delayed`` object.

        Parameters
        ----------
        optimize_graph : bool, optional
            If True [default], the graph is optimized before converting into
            ``dask.delayed`` objects.
        """
        from dask.delayed import Delayed
        dsk = self.__dask_graph__()
        if optimize_graph:
            dsk = self.__dask_optimize__(dsk, self.__dask_keys__())
        return Delayed(self.key, dsk)


class Bag(DaskMethodsMixin):
    """ Parallel collection of Python objects

    Examples
    --------
    Create Bag from sequence

    >>> import dask.bag as db
    >>> b = db.from_sequence(range(5))
    >>> list(b.filter(lambda x: x % 2 == 0).map(lambda x: x * 10))  # doctest: +SKIP
    [0, 20, 40]

    Create Bag from filename or globstring of filenames

    >>> b = db.read_text('/path/to/mydata.*.json.gz').map(json.loads)  # doctest: +SKIP

    Create manually (expert use)

    >>> dsk = {('x', 0): (range, 5),
    ...        ('x', 1): (range, 5),
    ...        ('x', 2): (range, 5)}
    >>> b = Bag(dsk, 'x', npartitions=3)

    >>> sorted(b.map(lambda x: x * 10))  # doctest: +SKIP
    [0, 0, 0, 10, 10, 10, 20, 20, 20, 30, 30, 30, 40, 40, 40]

    >>> int(b.fold(lambda x, y: x + y))  # doctest: +SKIP
    30
    """
    def __init__(self, dsk, name, npartitions):
        self.dask = dsk
        self.name = name
        self.npartitions = npartitions

    def __dask_graph__(self):
        return self.dask

    def __dask_keys__(self):
        return [(self.name, i) for i in range(self.npartitions)]

    def __dask_tokenize__(self):
        return self.name

    __dask_optimize__ = globalmethod(optimize, key='bag_optimize',
                                     falsey=dont_optimize)
    __dask_scheduler__ = staticmethod(mpget)

    def __dask_postcompute__(self):
        return finalize, ()

    def __dask_postpersist__(self):
        return type(self), (self.name, self.npartitions)

    def __str__(self):
        name = self.name if len(self.name) < 10 else self.name[:7] + '...'
        return 'dask.bag<%s, npartitions=%d>' % (name, self.npartitions)

    __repr__ = __str__

    str = property(fget=StringAccessor)

    def map(self, func, *args, **kwargs):
        """Apply a function elementwise across one or more bags.

        Note that all ``Bag`` arguments must be partitioned identically.

        Parameters
        ----------
        func : callable
        *args, **kwargs : Bag, Item, or object
            Extra arguments and keyword arguments to pass to ``func`` *after*
            the calling bag instance. Non-Bag args/kwargs are broadcasted
            across all calls to ``func``.

        Notes
        -----
        For calls with multiple `Bag` arguments, corresponding partitions
        should have the same length; if they do not, the call will error at
        compute time.

        Examples
        --------
        >>> import dask.bag as db
        >>> b = db.from_sequence(range(5), npartitions=2)
        >>> b2 = db.from_sequence(range(5, 10), npartitions=2)

        Apply a function to all elements in a bag:

        >>> b.map(lambda x: x + 1).compute()
        [1, 2, 3, 4, 5]

        Apply a function with arguments from multiple bags:

        >>> from operator import add
        >>> b.map(add, b2).compute()
        [5, 7, 9, 11, 13]

        Non-bag arguments are broadcast across all calls to the mapped
        function:

        >>> b.map(add, 1).compute()
        [1, 2, 3, 4, 5]

        Keyword arguments are also supported, and have the same semantics as
        regular arguments:

        >>> def myadd(x, y=0):
        ...     return x + y
        >>> b.map(myadd, y=b2).compute()
        [5, 7, 9, 11, 13]
        >>> b.map(myadd, y=1).compute()
        [1, 2, 3, 4, 5]

        Both arguments and keyword arguments can also be instances of
        ``dask.bag.Item``. Here we'll add the max value in the bag to each
        element:

        >>> b.map(myadd, b.max()).compute()
        [4, 5, 6, 7, 8]
        """
        return bag_map(func, self, *args, **kwargs)

    def starmap(self, func, **kwargs):
        """Apply a function using argument tuples from the given bag.

        This is similar to ``itertools.starmap``, except it also accepts
        keyword arguments. In pseudocode, this is could be written as:

        >>> def starmap(func, bag, **kwargs):
        ...     return (func(*args, **kwargs) for args in bag)

        Parameters
        ----------
        func : callable
        **kwargs : Item, Delayed, or object, optional
            Extra keyword arguments to pass to ``func``. These can either be
            normal objects, ``dask.bag.Item``, or ``dask.delayed.Delayed``.

        Examples
        --------
        >>> import dask.bag as db
        >>> data = [(1, 2), (3, 4), (5, 6), (7, 8), (9, 10)]
        >>> b = db.from_sequence(data, npartitions=2)

        Apply a function to each argument tuple:

        >>> from operator import add
        >>> b.starmap(add).compute()
        [3, 7, 11, 15, 19]

        Apply a function to each argument tuple, with additional keyword
        arguments:

        >>> def myadd(x, y, z=0):
        ...     return x + y + z
        >>> b.starmap(myadd, z=10).compute()
        [13, 17, 21, 25, 29]

        Keyword arguments can also be instances of ``dask.bag.Item`` or
        ``dask.delayed.Delayed``:

        >>> max_second = b.pluck(1).max()
        >>> max_second.compute()
        10
        >>> b.starmap(myadd, z=max_second).compute()
        [13, 17, 21, 25, 29]
        """
        name = 'starmap-{0}-{1}'.format(funcname(func),
                                        tokenize(self, func, kwargs))
        dsk = self.dask.copy()
        if kwargs:
            kw_dsk, kwargs = unpack_scalar_dask_kwargs(kwargs)
            dsk.update(kw_dsk)
        dsk.update({(name, i): (reify, (starmap_chunk, func, (self.name, i), kwargs))
                   for i in range(self.npartitions)})
        return type(self)(dsk, name, self.npartitions)

    @property
    def _args(self):
        return (self.dask, self.name, self.npartitions)

    def __getstate__(self):
        return self._args

    def __setstate__(self, state):
        self.dask, self.name, self.npartitions = state

    def filter(self, predicate):
        """ Filter elements in collection by a predicate function.

        >>> def iseven(x):
        ...     return x % 2 == 0

        >>> import dask.bag as db
        >>> b = db.from_sequence(range(5))
        >>> list(b.filter(iseven))  # doctest: +SKIP
        [0, 2, 4]
        """
        name = 'filter-{0}-{1}'.format(funcname(predicate),
                                       tokenize(self, predicate))
        dsk = dict(((name, i), (reify, (filter, predicate, (self.name, i))))
                   for i in range(self.npartitions))
        return type(self)(merge(self.dask, dsk), name, self.npartitions)

    def random_sample(self, prob, random_state=None):
        """ Return elements from bag with probability of ``prob``.

        Parameters
        ----------
        prob : float
            A float between 0 and 1, representing the probability that each
            element will be returned.
        random_state : int or random.Random, optional
            If an integer, will be used to seed a new ``random.Random`` object.
            If provided, results in deterministic sampling.

        Examples
        --------
        >>> import dask.bag as db
        >>> b = db.from_sequence(range(5))
        >>> list(b.random_sample(0.5, 42))
        [1, 4]
        >>> list(b.random_sample(0.5, 42))
        [1, 4]
        """
        if not 0 <= prob <= 1:
            raise ValueError('prob must be a number in the interval [0, 1]')
        if not isinstance(random_state, Random):
            random_state = Random(random_state)

        name = 'random-sample-%s' % tokenize(self, prob, random_state.getstate())
        state_data = random_state_data_python(self.npartitions, random_state)
        dsk = {(name, i): (reify, (random_sample, (self.name, i), state, prob))
               for i, state in zip(range(self.npartitions), state_data)}
        return type(self)(merge(self.dask, dsk), name, self.npartitions)

    def remove(self, predicate):
        """ Remove elements in collection that match predicate.

        >>> def iseven(x):
        ...     return x % 2 == 0

        >>> import dask.bag as db
        >>> b = db.from_sequence(range(5))
        >>> list(b.remove(iseven))  # doctest: +SKIP
        [1, 3]
        """
        name = 'remove-{0}-{1}'.format(funcname(predicate),
                                       tokenize(self, predicate))
        dsk = dict(((name, i), (reify, (remove, predicate, (self.name, i))))
                   for i in range(self.npartitions))
        return type(self)(merge(self.dask, dsk), name, self.npartitions)

    def map_partitions(self, func, *args, **kwargs):
        """Apply a function to every partition across one or more bags.

        Note that all ``Bag`` arguments must be partitioned identically.

        Parameters
        ----------
        func : callable
            The function to be called on every partition.
            This function should expect an ``Iterator`` or ``Iterable`` for
            every partition and should return an ``Iterator`` or ``Iterable``
            in return.
        *args, **kwargs : Bag, Item, Delayed, or object
            Arguments and keyword arguments to pass to ``func``.
            Partitions from this bag will be the first argument, and these will
            be passed *after*.

        Examples
        --------
        >>> import dask.bag as db
        >>> b = db.from_sequence(range(1, 101), npartitions=10)
        >>> def div(nums, den=1):
        ...     return [num / den for num in nums]

        Using a python object:

        >>> hi = b.max().compute()
        >>> hi
        100
        >>> b.map_partitions(div, den=hi).take(5)
        (0.01, 0.02, 0.03, 0.04, 0.05)

        Using an ``Item``:

        >>> b.map_partitions(div, den=b.max()).take(5)
        (0.01, 0.02, 0.03, 0.04, 0.05)

        Note that while both versions give the same output, the second forms a
        single graph, and then computes everything at once, and in some cases
        may be more efficient.
        """
        return map_partitions(func, self, *args, **kwargs)

    def pluck(self, key, default=no_default):
        """ Select item from all tuples/dicts in collection.

        >>> b = from_sequence([{'name': 'Alice', 'credits': [1, 2, 3]},
        ...                    {'name': 'Bob',   'credits': [10, 20]}])
        >>> list(b.pluck('name'))  # doctest: +SKIP
        ['Alice', 'Bob']
        >>> list(b.pluck('credits').pluck(0))  # doctest: +SKIP
        [1, 10]
        """
        name = 'pluck-' + tokenize(self, key, default)
        key = quote(key)
        if default == no_default:
            dsk = dict(((name, i), (list, (pluck, key, (self.name, i))))
                       for i in range(self.npartitions))
        else:
            dsk = dict(((name, i), (list, (pluck, key, (self.name, i), default)))
                       for i in range(self.npartitions))
        return type(self)(merge(self.dask, dsk), name, self.npartitions)

    def unzip(self, n):
        """Transform a bag of tuples to ``n`` bags of their elements.

        Examples
        --------
        >>> b = from_sequence([(i, i + 1, i + 2) for i in range(10)])
        >>> first, second, third = b.unzip(3)
        >>> isinstance(first, Bag)
        True
        >>> first.compute()
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        Note that this is equivalent to:

        >>> first, second, third = (b.pluck(i) for i in range(3))
        """
        return tuple(self.pluck(i) for i in range(n))

    @wraps(to_textfiles)
    def to_textfiles(self, path, name_function=None, compression='infer',
                     encoding=system_encoding, compute=True, get=None,
                     storage_options=None):
        return to_textfiles(self, path, name_function, compression, encoding,
                            compute, get=get, storage_options=storage_options)

    def fold(self, binop, combine=None, initial=no_default, split_every=None):
        """ Parallelizable reduction

        Fold is like the builtin function ``reduce`` except that it works in
        parallel.  Fold takes two binary operator functions, one to reduce each
        partition of our dataset and another to combine results between
        partitions

        1.  ``binop``: Binary operator to reduce within each partition
        2.  ``combine``:  Binary operator to combine results from binop

        Sequentially this would look like the following:

        >>> intermediates = [reduce(binop, part) for part in partitions]  # doctest: +SKIP
        >>> final = reduce(combine, intermediates)  # doctest: +SKIP

        If only one function is given then it is used for both functions
        ``binop`` and ``combine`` as in the following example to compute the
        sum:

        >>> def add(x, y):
        ...     return x + y

        >>> b = from_sequence(range(5))
        >>> b.fold(add).compute()  # doctest: +SKIP
        10

        In full form we provide both binary operators as well as their default
        arguments

        >>> b.fold(binop=add, combine=add, initial=0).compute()  # doctest: +SKIP
        10

        More complex binary operators are also doable

        >>> def add_to_set(acc, x):
        ...     ''' Add new element x to set acc '''
        ...     return acc | set([x])
        >>> b.fold(add_to_set, set.union, initial=set()).compute()  # doctest: +SKIP
        {1, 2, 3, 4, 5}

        See Also
        --------

        Bag.foldby
        """
        combine = combine or binop
        if initial is not no_default:
            return self.reduction(curry(_reduce, binop, initial=initial),
                                  curry(_reduce, combine),
                                  split_every=split_every)
        else:
            from toolz.curried import reduce
            return self.reduction(reduce(binop), reduce(combine),
                                  split_every=split_every)

    def frequencies(self, split_every=None):
        """ Count number of occurrences of each distinct element.

        >>> b = from_sequence(['Alice', 'Bob', 'Alice'])
        >>> dict(b.frequencies())  # doctest: +SKIP
        {'Alice': 2, 'Bob', 1}
        """
        return self.reduction(frequencies, merge_frequencies,
                              out_type=Bag, split_every=split_every,
                              name='frequencies').map_partitions(dictitems)

    def topk(self, k, key=None, split_every=None):
        """ K largest elements in collection

        Optionally ordered by some key function

        >>> b = from_sequence([10, 3, 5, 7, 11, 4])
        >>> list(b.topk(2))  # doctest: +SKIP
        [11, 10]

        >>> list(b.topk(2, lambda x: -x))  # doctest: +SKIP
        [3, 4]
        """
        if key:
            if callable(key) and takes_multiple_arguments(key):
                key = partial(apply, key)
            func = partial(topk, k, key=key)
        else:
            func = partial(topk, k)
        return self.reduction(func, compose(func, toolz.concat), out_type=Bag,
                              split_every=split_every, name='topk')

    def distinct(self):
        """ Distinct elements of collection

        Unordered without repeats.

        >>> b = from_sequence(['Alice', 'Bob', 'Alice'])
        >>> sorted(b.distinct())
        ['Alice', 'Bob']
        """
        return self.reduction(set, merge_distinct, out_type=Bag,
                              name='distinct')

    def reduction(self, perpartition, aggregate, split_every=None,
                  out_type=Item, name=None):
        """ Reduce collection with reduction operators.

        Parameters
        ----------
        perpartition: function
            reduction to apply to each partition
        aggregate: function
            reduction to apply to the results of all partitions
        split_every: int (optional)
            Group partitions into groups of this size while performing reduction
            Defaults to 8
        out_type: {Bag, Item}
            The out type of the result, Item if a single element, Bag if a list
            of elements.  Defaults to Item.

        Examples
        --------
        >>> b = from_sequence(range(10))
        >>> b.reduction(sum, sum).compute()
        45
        """
        if split_every is None:
            split_every = 8
        if split_every is False:
            split_every = self.npartitions

        token = tokenize(self, perpartition, aggregate, split_every)
        a = '%s-part-%s' % (name or funcname(perpartition), token)
        is_last = self.npartitions == 1
        dsk = {(a, i): (empty_safe_apply, perpartition, (self.name, i), is_last)
               for i in range(self.npartitions)}
        k = self.npartitions
        b = a
        fmt = '%s-aggregate-%s' % (name or funcname(aggregate), token)
        depth = 0

        while k > split_every:
            c = fmt + str(depth)
            dsk2 = dict(((c, i), (empty_safe_aggregate, aggregate,
                                  [(b, j) for j in inds], False))
                        for i, inds in enumerate(partition_all(split_every,
                                                               range(k))))
            dsk.update(dsk2)
            k = len(dsk2)
            b = c
            depth += 1

        dsk[(fmt, 0)] = (empty_safe_aggregate, aggregate,
                         [(b, j) for j in range(k)], True)

        if out_type is Item:
            dsk[fmt] = dsk.pop((fmt, 0))
            return Item(merge(self.dask, dsk), fmt)
        else:
            return Bag(merge(self.dask, dsk), fmt, 1)

    def sum(self, split_every=None):
        """ Sum all elements """
        return self.reduction(sum, sum, split_every=split_every)

    def max(self, split_every=None):
        """ Maximum element """
        return self.reduction(max, max, split_every=split_every)

    def min(self, split_every=None):
        """ Minimum element """
        return self.reduction(min, min, split_every=split_every)

    def any(self, split_every=None):
        """ Are any of the elements truthy? """
        return self.reduction(any, any, split_every=split_every)

    def all(self, split_every=None):
        """ Are all elements truthy? """
        return self.reduction(all, all, split_every=split_every)

    def count(self, split_every=None):
        """ Count the number of elements. """
        return self.reduction(count, sum, split_every=split_every)

    def mean(self):
        """ Arithmetic mean """
        def mean_chunk(seq):
            total, n = 0.0, 0
            for x in seq:
                total += x
                n += 1
            return total, n

        def mean_aggregate(x):
            totals, counts = list(zip(*x))
            return 1.0 * sum(totals) / sum(counts)

        return self.reduction(mean_chunk, mean_aggregate, split_every=False)

    def var(self, ddof=0):
        """ Variance """
        def var_chunk(seq):
            squares, total, n = 0.0, 0.0, 0
            for x in seq:
                squares += x**2
                total += x
                n += 1
            return squares, total, n

        def var_aggregate(x):
            squares, totals, counts = list(zip(*x))
            x2, x, n = float(sum(squares)), float(sum(totals)), sum(counts)
            result = (x2 / n) - (x / n)**2
            return result * n / (n - ddof)

        return self.reduction(var_chunk, var_aggregate, split_every=False)

    def std(self, ddof=0):
        """ Standard deviation """
        return self.var(ddof=ddof).apply(math.sqrt)

    def join(self, other, on_self, on_other=None):
        """ Joins collection with another collection.

        Other collection must be one of the following:

        1.  An iterable.  We recommend tuples over lists for internal
            performance reasons.
        2.  A delayed object, pointing to a tuple.  This is recommended if the
            other collection is sizable and you're using the distributed
            scheduler.  Dask is able to pass around data wrapped in delayed
            objects with greater sophistication.
        3.  A Bag with a single partition

        You might also consider Dask Dataframe, whose join operations are much
        more heavily optimized.

        Parameters
        ----------
        other: Iterable, Delayed, Bag
            Other collection on which to join
        on_self: callable
            Function to call on elements in this collection to determine a
            match
        on_other: callable (defaults to on_self)
            Function to call on elements in the other collection to determine a
            match

        Examples
        --------
        >>> people = from_sequence(['Alice', 'Bob', 'Charlie'])
        >>> fruit = ['Apple', 'Apricot', 'Banana']
        >>> list(people.join(fruit, lambda x: x[0]))  # doctest: +SKIP
        [('Apple', 'Alice'), ('Apricot', 'Alice'), ('Banana', 'Bob')]
        """
        name = 'join-' + tokenize(self, other, on_self, on_other)
        dsk = {}
        if isinstance(other, Bag):
            if other.npartitions == 1:
                dsk.update(other.dask)
                other = other.__dask_keys__()[0]
                dsk['join-%s-other' % name] = (list, other)
            else:
                msg = ("Multi-bag joins are not implemented. "
                       "We recommend Dask dataframe if appropriate")
                raise NotImplementedError(msg)
        elif isinstance(other, Delayed):
            dsk.update(other.dask)
            other = other._key
        elif isinstance(other, Iterable):
            other = other
        else:
            msg = ("Joined argument must be single-partition Bag, "
                   " delayed object, or Iterable, got %s" %
                   type(other).__name)
            raise TypeError(msg)

        if on_other is None:
            on_other = on_self

        dsk.update({(name, i): (list, (join, on_other, other,
                                       on_self, (self.name, i)))
                   for i in range(self.npartitions)})
        return type(self)(merge(self.dask, dsk), name, self.npartitions)

    def product(self, other):
        """ Cartesian product between two bags. """
        assert isinstance(other, Bag)
        name = 'product-' + tokenize(self, other)
        n, m = self.npartitions, other.npartitions
        dsk = dict(((name, i * m + j),
                   (list, (itertools.product, (self.name, i),
                                              (other.name, j))))
                   for i in range(n) for j in range(m))
        return type(self)(merge(self.dask, other.dask, dsk), name, n * m)

    def foldby(self, key, binop, initial=no_default, combine=None,
               combine_initial=no_default, split_every=None):
        """ Combined reduction and groupby.

        Foldby provides a combined groupby and reduce for efficient parallel
        split-apply-combine tasks.

        The computation

        >>> b.foldby(key, binop, init)                        # doctest: +SKIP

        is equivalent to the following:

        >>> def reduction(group):                               # doctest: +SKIP
        ...     return reduce(binop, group, init)               # doctest: +SKIP

        >>> b.groupby(key).map(lambda (k, v): (k, reduction(v)))# doctest: +SKIP

        But uses minimal communication and so is *much* faster.

        >>> b = from_sequence(range(10))
        >>> iseven = lambda x: x % 2 == 0
        >>> add = lambda x, y: x + y
        >>> dict(b.foldby(iseven, add))                         # doctest: +SKIP
        {True: 20, False: 25}

        **Key Function**

        The key function determines how to group the elements in your bag.
        In the common case where your bag holds dictionaries then the key
        function often gets out one of those elements.

        >>> def key(x):
        ...     return x['name']

        This case is so common that it is special cased, and if you provide a
        key that is not a callable function then dask.bag will turn it into one
        automatically.  The following are equivalent:

        >>> b.foldby(lambda x: x['name'], ...)  # doctest: +SKIP
        >>> b.foldby('name', ...)  # doctest: +SKIP

        **Binops**

        It can be tricky to construct the right binary operators to perform
        analytic queries.  The ``foldby`` method accepts two binary operators,
        ``binop`` and ``combine``.  Binary operators two inputs and output must
        have the same type.

        Binop takes a running total and a new element and produces a new total:

        >>> def binop(total, x):
        ...     return total + x['amount']

        Combine takes two totals and combines them:

        >>> def combine(total1, total2):
        ...     return total1 + total2

        Each of these binary operators may have a default first value for
        total, before any other value is seen.  For addition binary operators
        like above this is often ``0`` or the identity element for your
        operation.

        **split_every**

        Group partitions into groups of this size while performing reduction.
        Defaults to 8.

        >>> b.foldby('name', binop, 0, combine, 0)  # doctest: +SKIP

        See Also
        --------

        toolz.reduceby
        pyspark.combineByKey
        """
        if split_every is None:
            split_every = 8
        if split_every is False:
            split_every = self.npartitions

        token = tokenize(self, key, binop, initial, combine, combine_initial)
        a = 'foldby-a-' + token
        if combine is None:
            combine = binop
        if initial is not no_default:
            dsk = {(a, i): (reduceby, key, binop, (self.name, i), initial)
                   for i in range(self.npartitions)}
        else:
            dsk = {(a, i): (reduceby, key, binop, (self.name, i))
                   for i in range(self.npartitions)}

        def combine2(acc, x):
            return combine(acc, x[1])

        depth = 0
        k = self.npartitions
        b = a
        while k > split_every:
            c = b + str(depth)
            if combine_initial is not no_default:
                dsk2 = {(c, i): (reduceby, 0, combine2,
                                 (toolz.concat, (map, dictitems,
                                                 [(b, j) for j in inds])),
                                 combine_initial)
                        for i, inds in enumerate(partition_all(split_every,
                                                               range(k)))}
            else:
                dsk2 = {(c, i): (merge_with, (partial, reduce, combine),
                                 [(b, j) for j in inds])
                        for i, inds in enumerate(partition_all(split_every,
                                                               range(k)))}
            dsk.update(dsk2)
            k = len(dsk2)
            b = c
            depth += 1

        e = 'foldby-b-' + token
        if combine_initial is not no_default:
            dsk[(e, 0)] = (dictitems, (reduceby, 0, combine2,
                                       (toolz.concat, (map, dictitems,
                                                       [(b, j) for j in range(k)])),
                                       combine_initial))
        else:
            dsk[(e, 0)] = (dictitems, (merge_with, (partial, reduce, combine),
                                       [(b, j) for j in range(k)]))

        return type(self)(merge(self.dask, dsk), e, 1)

    def take(self, k, npartitions=1, compute=True, warn=True):
        """ Take the first k elements.

        Parameters
        ----------
        k : int
            The number of elements to return
        npartitions : int, optional
            Elements are only taken from the first ``npartitions``, with a
            default of 1. If there are fewer than ``k`` rows in the first
            ``npartitions`` a warning will be raised and any found rows
            returned. Pass -1 to use all partitions.
        compute : bool, optional
            Whether to compute the result, default is True.
        warn : bool, optional
            Whether to warn if the number of elements returned is less than
            requested, default is True.

        >>> b = from_sequence(range(10))
        >>> b.take(3)  # doctest: +SKIP
        (0, 1, 2)
        """

        if npartitions <= -1:
            npartitions = self.npartitions
        if npartitions > self.npartitions:
            raise ValueError("only {} partitions, take "
                             "received {}".format(self.npartitions, npartitions))

        token = tokenize(self, k, npartitions)
        name = 'take-' + token

        if npartitions > 1:
            name_p = 'take-partial-' + token

            dsk = {}
            for i in range(npartitions):
                dsk[(name_p, i)] = (list, (take, k, (self.name, i)))

            concat = (toolz.concat, ([(name_p, i) for i in range(npartitions)]))
            dsk[(name, 0)] = (safe_take, k, concat, warn)
        else:
            dsk = {(name, 0): (safe_take, k, (self.name, 0), warn)}

        b = Bag(merge(self.dask, dsk), name, 1)

        if compute:
            return tuple(b.compute())
        else:
            return b

    def flatten(self):
        """ Concatenate nested lists into one long list.

        >>> b = from_sequence([[1], [2, 3]])
        >>> list(b)
        [[1], [2, 3]]

        >>> list(b.flatten())
        [1, 2, 3]
        """
        name = 'flatten-' + tokenize(self)
        dsk = dict(((name, i), (list, (toolz.concat, (self.name, i))))
                   for i in range(self.npartitions))
        return type(self)(merge(self.dask, dsk), name, self.npartitions)

    def __iter__(self):
        return iter(self.compute())

    def groupby(self, grouper, method=None, npartitions=None, blocksize=2**20,
                max_branch=None):
        """ Group collection by key function

        This requires a full dataset read, serialization and shuffle.
        This is expensive.  If possible you should use ``foldby``.

        Parameters
        ----------
        grouper: function
            Function on which to group elements
        method: str
            Either 'disk' for an on-disk shuffle or 'tasks' to use the task
            scheduling framework.  Use 'disk' if you are on a single machine
            and 'tasks' if you are on a distributed cluster.
        npartitions: int
            If using the disk-based shuffle, the number of output partitions
        blocksize: int
            If using the disk-based shuffle, the size of shuffle blocks (bytes)
        max_branch: int
            If using the task-based shuffle, the amount of splitting each
            partition undergoes.  Increase this for fewer copies but more
            scheduler overhead.

        Examples
        --------
        >>> b = from_sequence(range(10))
        >>> iseven = lambda x: x % 2 == 0
        >>> dict(b.groupby(iseven))  # doctest: +SKIP
        {True: [0, 2, 4, 6, 8], False: [1, 3, 5, 7, 9]}

        See Also
        --------
        Bag.foldby
        """
        if method is None:
            get = _globals.get('get')
            if (isinstance(get, types.MethodType) and
               'distributed' in get.__func__.__module__):
                method = 'tasks'
            else:
                method = 'disk'
        if method == 'disk':
            return groupby_disk(self, grouper, npartitions=npartitions,
                                blocksize=blocksize)
        elif method == 'tasks':
            return groupby_tasks(self, grouper, max_branch=max_branch)
        else:
            msg = "Shuffle method must be 'disk' or 'tasks'"
            raise NotImplementedError(msg)

    def to_dataframe(self, meta=None, columns=None):
        """ Create Dask Dataframe from a Dask Bag.

        Bag should contain tuples, dict records, or scalars.

        Index will not be particularly meaningful.  Use ``reindex`` afterwards
        if necessary.

        Parameters
        ----------
        meta : pd.DataFrame, dict, iterable, optional
            An empty ``pd.DataFrame`` that matches the dtypes and column names
            of the output. This metadata is necessary for many algorithms in
            dask dataframe to work.  For ease of use, some alternative inputs
            are also available. Instead of a ``DataFrame``, a ``dict`` of
            ``{name: dtype}`` or iterable of ``(name, dtype)`` can be provided.
            If not provided or a list, a single element from the first
            partition will be computed, triggering a potentially expensive call
            to ``compute``. This may lead to unexpected results, so providing
            ``meta`` is recommended. For more information, see
            ``dask.dataframe.utils.make_meta``.
        columns : sequence, optional
            Column names to use. If the passed data do not have names
            associated with them, this argument provides names for the columns.
            Otherwise this argument indicates the order of the columns in the
            result (any names not found in the data will become all-NA
            columns).  Note that if ``meta`` is provided, column names will be
            taken from there and this parameter is invalid.

        Examples
        --------
        >>> import dask.bag as db
        >>> b = db.from_sequence([{'name': 'Alice',   'balance': 100},
        ...                       {'name': 'Bob',     'balance': 200},
        ...                       {'name': 'Charlie', 'balance': 300}],
        ...                      npartitions=2)
        >>> df = b.to_dataframe()

        >>> df.compute()
           balance     name
        0      100    Alice
        1      200      Bob
        0      300  Charlie
        """
        import pandas as pd
        import dask.dataframe as dd
        if meta is None:
            if isinstance(columns, pd.DataFrame):
                warnings.warn("Passing metadata to `columns` is deprecated. "
                              "Please use the `meta` keyword instead.")
                meta = columns
            else:
                head = self.take(1, warn=False)
                if len(head) == 0:
                    raise ValueError("`dask.bag.Bag.to_dataframe` failed to "
                                     "properly infer metadata, please pass in "
                                     "metadata via the `meta` keyword")
                meta = pd.DataFrame(list(head), columns=columns)
        elif columns is not None:
            raise ValueError("Can't specify both `meta` and `columns`")
        else:
            meta = dd.utils.make_meta(meta)
        # Serializing the columns and dtypes is much smaller than serializing
        # the empty frame
        cols = list(meta.columns)
        dtypes = meta.dtypes.to_dict()
        name = 'to_dataframe-' + tokenize(self, cols, dtypes)
        dsk = self.__dask_optimize__(self.dask, self.__dask_keys__())
        dsk.update({(name, i): (to_dataframe, (self.name, i), cols, dtypes)
                    for i in range(self.npartitions)})
        divisions = [None] * (self.npartitions + 1)
        return dd.DataFrame(dsk, name, meta, divisions)

    def to_delayed(self, optimize_graph=True):
        """Convert into a list of ``dask.delayed`` objects, one per partition.

        Parameters
        ----------
        optimize_graph : bool, optional
            If True [default], the graph is optimized before converting into
            ``dask.delayed`` objects.

        See Also
        --------
        dask.bag.from_delayed
        """
        from dask.delayed import Delayed
        keys = self.__dask_keys__()
        dsk = self.__dask_graph__()
        if optimize_graph:
            dsk = self.__dask_optimize__(dsk, keys)
        return [Delayed(k, dsk) for k in keys]

    def repartition(self, npartitions):
        """ Coalesce bag into fewer partitions.

        Examples
        --------
        >>> b.repartition(5)  # set to have 5 partitions  # doctest: +SKIP
        """
        new_name = 'repartition-%d-%s' % (npartitions, tokenize(self, npartitions))
        if npartitions == self.npartitions:
            return self
        elif npartitions < self.npartitions:
            ratio = self.npartitions / npartitions
            new_partitions_boundaries = [int(old_partition_index * ratio)
                                         for old_partition_index in range(npartitions + 1)]

            dsk = {}
            for new_partition_index in range(npartitions):
                value = (list, (toolz.concat,
                                [(self.name, old_partition_index)
                                 for old_partition_index in
                                 range(new_partitions_boundaries[new_partition_index],
                                       new_partitions_boundaries[new_partition_index + 1])]))
                dsk[new_name, new_partition_index] = value
        else:  # npartitions > self.npartitions
            ratio = npartitions / self.npartitions
            split_name = 'split-%s' % tokenize(self, npartitions)
            dsk = {}
            last = 0
            j = 0
            for i in range(self.npartitions):
                new = last + ratio
                if i == self.npartitions - 1:
                    k = npartitions - j
                else:
                    k = int(new - last)
                dsk[(split_name, i)] = (split, (self.name, i), k)
                for jj in range(k):
                    dsk[(new_name, j)] = (getitem, (split_name, i), jj)
                    j += 1
                last = new

        return Bag(dsk=merge(self.dask, dsk), name=new_name, npartitions=npartitions)

    def accumulate(self, binop, initial=no_default):
        """ Repeatedly apply binary function to a sequence, accumulating results.

        This assumes that the bag is ordered.  While this is typically the case
        not all Dask.bag functions preserve this property.

        Examples
        --------
        >>> from operator import add
        >>> b = from_sequence([1, 2, 3, 4, 5], npartitions=2)
        >>> b.accumulate(add).compute()  # doctest: +SKIP
        [1, 3, 6, 10, 15]

        Accumulate also takes an optional argument that will be used as the
        first value.

        >>> b.accumulate(add, initial=-1)  # doctest: +SKIP
        [-1, 0, 2, 5, 9, 14]
        """
        if not _implement_accumulate:
            raise NotImplementedError("accumulate requires `toolz` > 0.7.4"
                                      " or `cytoolz` > 0.7.3.")
        token = tokenize(self, binop, initial)
        binop_name = funcname(binop)
        a = '%s-part-%s' % (binop_name, token)
        b = '%s-first-%s' % (binop_name, token)
        c = '%s-second-%s' % (binop_name, token)
        dsk = {(a, 0): (accumulate_part, binop, (self.name, 0), initial, True),
               (b, 0): (first, (a, 0)),
               (c, 0): (second, (a, 0))}
        for i in range(1, self.npartitions):
            dsk[(a, i)] = (accumulate_part, binop, (self.name, i), (c, i - 1))
            dsk[(b, i)] = (first, (a, i))
            dsk[(c, i)] = (second, (a, i))
        return Bag(merge(self.dask, dsk), b, self.npartitions)


def accumulate_part(binop, seq, initial, is_first=False):
    if initial == no_default:
        res = list(accumulate(binop, seq))
    else:
        res = list(accumulate(binop, seq, initial=initial))
    if is_first:
        return res, res[-1] if res else [], initial
    return res[1:], res[-1]


def partition(grouper, sequence, npartitions, p, nelements=2**20):
    """ Partition a bag along a grouper, store partitions on disk. """
    for block in partition_all(nelements, sequence):
        d = groupby(grouper, block)
        d2 = defaultdict(list)
        for k, v in d.items():
            d2[abs(hash(k)) % npartitions].extend(v)
        p.append(d2, fsync=True)
    return p


def collect(grouper, group, p, barrier_token):
    """ Collect partitions from disk and yield k,v group pairs. """
    d = groupby(grouper, p.get(group, lock=False))
    return list(d.items())


def from_sequence(seq, partition_size=None, npartitions=None):
    """ Create a dask Bag from Python sequence.

    This sequence should be relatively small in memory.  Dask Bag works
    best when it handles loading your data itself.  Commonly we load a
    sequence of filenames into a Bag and then use ``.map`` to open them.

    Parameters
    ----------
    seq: Iterable
        A sequence of elements to put into the dask
    partition_size: int (optional)
        The length of each partition
    npartitions: int (optional)
        The number of desired partitions

    It is best to provide either ``partition_size`` or ``npartitions``
    (though not both.)

    Examples
    --------
    >>> b = from_sequence(['Alice', 'Bob', 'Chuck'], partition_size=2)

    See Also
    --------
    read_text: Create bag from text files
    """
    seq = list(seq)
    if npartitions and not partition_size:
        partition_size = int(math.ceil(len(seq) / npartitions))
    if npartitions is None and partition_size is None:
        if len(seq) < 100:
            partition_size = 1
        else:
            partition_size = int(len(seq) / 100)

    parts = list(partition_all(partition_size, seq))
    name = 'from_sequence-' + tokenize(seq, partition_size)
    d = dict(((name, i), list(part)) for i, part in enumerate(parts))
    return Bag(d, name, len(d))


def from_url(urls):
    """Create a dask Bag from a url.

    Examples
    --------
    >>> a = from_url('http://raw.githubusercontent.com/dask/dask/master/README.rst')  # doctest: +SKIP
    >>> a.npartitions  # doctest: +SKIP
    1

    >>> a.take(8)  # doctest: +SKIP
    (b'Dask\\n',
     b'====\\n',
     b'\\n',
     b'|Build Status| |Coverage| |Doc Status| |Gitter| |Version Status|\\n',
     b'\\n',
     b'Dask is a flexible parallel computing library for analytics.  See\\n',
     b'documentation_ for more information.\\n',
     b'\\n')

    >>> b = from_url(['http://github.com', 'http://google.com'])  # doctest: +SKIP
    >>> b.npartitions  # doctest: +SKIP
    2
    """
    if isinstance(urls, str):
        urls = [urls]
    name = 'from_url-' + uuid.uuid4().hex
    dsk = {}
    for i, u in enumerate(urls):
        dsk[(name, i)] = (list, (urlopen, u))
    return Bag(dsk, name, len(urls))


def dictitems(d):
    """ A pickleable version of dict.items

    >>> dictitems({'x': 1})
    [('x', 1)]
    """
    return list(d.items())


def concat(bags):
    """ Concatenate many bags together, unioning all elements.

    >>> import dask.bag as db
    >>> a = db.from_sequence([1, 2, 3])
    >>> b = db.from_sequence([4, 5, 6])
    >>> c = db.concat([a, b])

    >>> list(c)
    [1, 2, 3, 4, 5, 6]
    """
    name = 'concat-' + tokenize(*bags)
    counter = itertools.count(0)
    dsk = {(name, next(counter)): key
           for bag in bags for key in bag.__dask_keys__()}
    return Bag(merge(dsk, *[b.dask for b in bags]), name, len(dsk))


def reify(seq):
    if isinstance(seq, Iterator):
        seq = list(seq)
    if seq and isinstance(seq[0], Iterator):
        seq = list(map(list, seq))
    return seq


def from_delayed(values):
    """ Create bag from many dask Delayed objects.

    These objects will become the partitions of the resulting Bag.  They should
    evaluate to a ``list`` or some other concrete sequence.

    Parameters
    ----------
    values: list of delayed values
        An iterable of dask Delayed objects.  Each evaluating to a list.

    Returns
    -------
    Bag

    Examples
    --------
    >>> x, y, z = [delayed(load_sequence_from_file)(fn)
    ...             for fn in filenames] # doctest: +SKIP
    >>> b = from_delayed([x, y, z])  # doctest: +SKIP

    See also
    --------
    dask.delayed
    """
    from dask.delayed import Delayed, delayed
    if isinstance(values, Delayed):
        values = [values]
    values = [delayed(v)
              if not isinstance(v, Delayed) and hasattr(v, 'key')
              else v
              for v in values]
    dsk = merge(ensure_dict(v.dask) for v in values)

    name = 'bag-from-delayed-' + tokenize(*values)
    names = [(name, i) for i in range(len(values))]
    values = [(reify, v.key) for v in values]
    dsk2 = dict(zip(names, values))

    return Bag(merge(dsk, dsk2), name, len(values))


def merge_distinct(seqs):
    return set().union(*seqs)


def merge_frequencies(seqs):
    if isinstance(seqs, Iterable):
        seqs = list(seqs)
    if not seqs:
        return {}
    first, rest = seqs[0], seqs[1:]
    if not rest:
        return first
    out = defaultdict(int)
    out.update(first)
    for d in rest:
        for k, v in iteritems(d):
            out[k] += v
    return out


def bag_range(n, npartitions):
    """ Numbers from zero to n

    Examples
    --------

    >>> import dask.bag as db
    >>> b = db.range(5, npartitions=2)
    >>> list(b)
    [0, 1, 2, 3, 4]
    """
    size = n // npartitions
    name = 'range-%d-npartitions-%d' % (n, npartitions)
    ijs = list(enumerate(take(npartitions, range(0, n, size))))
    dsk = dict(((name, i), (reify, (range, j, min(j + size, n))))
               for i, j in ijs)

    if n % npartitions != 0:
        i, j = ijs[-1]
        dsk[(name, i)] = (reify, (range, j, n))

    return Bag(dsk, name, npartitions)


def bag_zip(*bags):
    """ Partition-wise bag zip

    All passed bags must have the same number of partitions.

    NOTE: corresponding partitions should have the same length; if they do not,
    the "extra" elements from the longer partition(s) will be dropped.  If you
    have this case chances are that what you really need is a data alignment
    mechanism like pandas's, and not a missing value filler like zip_longest.

    Examples
    --------

    Correct usage:

    >>> import dask.bag as db
    >>> evens = db.from_sequence(range(0, 10, 2), partition_size=4)
    >>> odds = db.from_sequence(range(1, 10, 2), partition_size=4)
    >>> pairs = db.zip(evens, odds)
    >>> list(pairs)
    [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9)]

    Incorrect usage:

    >>> numbers = db.range(20) # doctest: +SKIP
    >>> fizz = numbers.filter(lambda n: n % 3 == 0) # doctest: +SKIP
    >>> buzz = numbers.filter(lambda n: n % 5 == 0) # doctest: +SKIP
    >>> fizzbuzz = db.zip(fizz, buzz) # doctest: +SKIP
    >>> list(fizzbuzzz) # doctest: +SKIP
    [(0, 0), (3, 5), (6, 10), (9, 15), (12, 20), (15, 25), (18, 30)]

    When what you really wanted was more along the lines of the following:

    >>> list(fizzbuzzz) # doctest: +SKIP
    [(0, 0), (3, None), (None, 5), (6, None), (None 10), (9, None),
    (12, None), (15, 15), (18, None), (None, 20), (None, 25), (None, 30)]
    """
    npartitions = bags[0].npartitions
    assert all(bag.npartitions == npartitions for bag in bags)
    # TODO: do more checks

    name = 'zip-' + tokenize(*bags)
    dsk = dict(
        ((name, i), (reify, (zip,) + tuple((bag.name, i) for bag in bags)))
        for i in range(npartitions))
    bags_dsk = merge(*(bag.dask for bag in bags))
    return Bag(merge(bags_dsk, dsk), name, npartitions)


def map_chunk(f, args, bag_kwargs, kwargs):
    if kwargs:
        f = partial(f, **kwargs)

    args = [iter(a) for a in args]
    iters = list(args)
    if bag_kwargs:
        keys = list(bag_kwargs)
        kw_val_iters = [iter(v) for v in bag_kwargs.values()]
        iters.extend(kw_val_iters)
        kw_iter = (dict(zip(keys, k)) for k in zip(*kw_val_iters))
        if args:
            for a, k in zip(zip(*args), kw_iter):
                yield f(*a, **k)
        else:
            for k in kw_iter:
                yield f(**k)
    else:
        for a in zip(*args):
            yield f(*a)

    # Check that all iterators are fully exhausted
    if len(iters) > 1:
        for i in iters:
            if isinstance(i, itertools.repeat):
                continue
            try:
                next(i)
            except StopIteration:
                pass
            else:
                msg = ("map called with multiple bags that aren't identically "
                       "partitioned. Please ensure that all bag arguments "
                       "have the same partition lengths")
                raise ValueError(msg)


def starmap_chunk(f, x, kwargs):
    if kwargs:
        f = partial(f, **kwargs)
    return itertools.starmap(f, x)


def unpack_scalar_dask_kwargs(kwargs):
    """Extracts dask values from kwargs.

    Currently only ``dask.bag.Item`` and ``dask.delayed.Delayed`` are
    supported.  Returns a merged dask graph and a task resulting in a keyword
    dict.
    """
    dsk = {}
    kwargs2 = {}
    for k, v in kwargs.items():
        if isinstance(v, (Delayed, Item)):
            dsk.update(ensure_dict(v.dask))
            kwargs2[k] = v.key
        elif is_dask_collection(v):
            raise NotImplementedError("dask.bag doesn't support kwargs of "
                                      "type %s" % type(v).__name__)
        else:
            kwargs2[k] = v
    if dsk:
        kwargs = (dict, (zip, list(kwargs2), list(kwargs2.values())))
    return dsk, kwargs


def bag_map(func, *args, **kwargs):
    """Apply a function elementwise across one or more bags.

    Note that all ``Bag`` arguments must be partitioned identically.

    Parameters
    ----------
    func : callable
    *args, **kwargs : Bag, Item, Delayed, or object
        Arguments and keyword arguments to pass to ``func``. Non-Bag args/kwargs
        are broadcasted across all calls to ``func``.

    Notes
    -----
    For calls with multiple `Bag` arguments, corresponding partitions should
    have the same length; if they do not, the call will error at compute time.

    Examples
    --------
    >>> import dask.bag as db
    >>> b = db.from_sequence(range(5), npartitions=2)
    >>> b2 = db.from_sequence(range(5, 10), npartitions=2)

    Apply a function to all elements in a bag:

    >>> db.map(lambda x: x + 1, b).compute()
    [1, 2, 3, 4, 5]

    Apply a function with arguments from multiple bags:

    >>> from operator import add
    >>> db.map(add, b, b2).compute()
    [5, 7, 9, 11, 13]

    Non-bag arguments are broadcast across all calls to the mapped function:

    >>> db.map(add, b, 1).compute()
    [1, 2, 3, 4, 5]

    Keyword arguments are also supported, and have the same semantics as
    regular arguments:

    >>> def myadd(x, y=0):
    ...     return x + y
    >>> db.map(myadd, b, y=b2).compute()
    [5, 7, 9, 11, 13]
    >>> db.map(myadd, b, y=1).compute()
    [1, 2, 3, 4, 5]

    Both arguments and keyword arguments can also be instances of
    ``dask.bag.Item`` or ``dask.delayed.Delayed``. Here we'll add the max value
    in the bag to each element:

    >>> db.map(myadd, b, b.max()).compute()
    [4, 5, 6, 7, 8]
    """
    name = 'map-%s-%s' % (funcname(func), tokenize(func, args, kwargs))
    dsk = {}

    bags = []
    args2 = []
    for a in args:
        if isinstance(a, Bag):
            bags.append(a)
            args2.append(a)
            dsk.update(a.dask)
        elif isinstance(a, (Item, Delayed)):
            args2.append((itertools.repeat, a.key))
            dsk.update(ensure_dict(a.dask))
        else:
            args2.append((itertools.repeat, a))

    bag_kwargs = {}
    other_kwargs = {}
    for k, v in kwargs.items():
        if isinstance(v, Bag):
            bag_kwargs[k] = v
            bags.append(v)
            dsk.update(v.dask)
        else:
            other_kwargs[k] = v

    kw_dsk, other_kwargs = unpack_scalar_dask_kwargs(other_kwargs)
    dsk.update(kw_dsk)

    if not bags:
        raise ValueError("At least one argument must be a Bag.")

    npartitions = {b.npartitions for b in bags}
    if len(npartitions) > 1:
        raise ValueError("All bags must have the same number of partitions.")
    npartitions = npartitions.pop()

    def build_args(n):
        return [(a.name, n) if isinstance(a, Bag) else a for a in args2]

    def build_bag_kwargs(n):
        if not bag_kwargs:
            return None
        return (dict, (zip, list(bag_kwargs),
                       [(b.name, n) for b in bag_kwargs.values()]))

    dsk.update({(name, n): (reify, (map_chunk, func, build_args(n),
                                    build_bag_kwargs(n), other_kwargs))
                for n in range(npartitions)})

    # If all bags are the same type, use that type, otherwise fallback to Bag
    return_type = set(map(type, bags))
    return_type = return_type.pop() if len(return_type) == 1 else Bag

    return return_type(dsk, name, npartitions)


def map_partitions(func, *args, **kwargs):
    """Apply a function to every partition across one or more bags.

    Note that all ``Bag`` arguments must be partitioned identically.

    Parameters
    ----------
    func : callable
    *args, **kwargs : Bag, Item, Delayed, or object
        Arguments and keyword arguments to pass to ``func``.

    Examples
    --------
    >>> import dask.bag as db
    >>> b = db.from_sequence(range(1, 101), npartitions=10)
    >>> def div(nums, den=1):
    ...     return [num / den for num in nums]

    Using a python object:

    >>> hi = b.max().compute()
    >>> hi
    100
    >>> b.map_partitions(div, den=hi).take(5)
    (0.01, 0.02, 0.03, 0.04, 0.05)

    Using an ``Item``:

    >>> b.map_partitions(div, den=b.max()).take(5)
    (0.01, 0.02, 0.03, 0.04, 0.05)

    Note that while both versions give the same output, the second forms a
    single graph, and then computes everything at once, and in some cases
    may be more efficient.
    """
    name = 'map-partitions-%s-%s' % (funcname(func),
                                     tokenize(func, args, kwargs))

    # Extract bag arguments, build initial graph
    bags = []
    dsk = {}
    for vals in [args, kwargs.values()]:
        for a in vals:
            if isinstance(a, (Bag, Item, Delayed)):
                dsk.update(ensure_dict(a.dask))
                if isinstance(a, Bag):
                    bags.append(a)
            elif is_dask_collection(a):
                raise NotImplementedError("dask.bag doesn't support args of "
                                          "type %s" % type(a).__name__)

    if not bags:
        raise ValueError("At least one argument must be a Bag.")

    npartitions = {b.npartitions for b in bags}
    if len(npartitions) > 1:
        raise ValueError("All bags must have the same number of partitions.")
    npartitions = npartitions.pop()

    def build_task(n):
        args2 = [(a.name, n) if isinstance(a, Bag) else a.key
                 if isinstance(a, (Item, Delayed)) else a for a in args]

        if any(isinstance(v, (Bag, Item, Delayed)) for v in kwargs.values()):
            vals = [(v.name, n) if isinstance(v, Bag) else v.key
                    if isinstance(v, (Item, Delayed)) else v
                    for v in kwargs.values()]
            kwargs2 = (dict, (zip, list(kwargs), vals))
        else:
            kwargs2 = kwargs

        if kwargs2 or len(args2) > 1:
            return (apply, func, args2, kwargs2)
        return (func, args2[0])

    dsk.update({(name, n): build_task(n) for n in range(npartitions)})

    # If all bags are the same type, use that type, otherwise fallback to Bag
    return_type = set(map(type, bags))
    return_type = return_type.pop() if len(return_type) == 1 else Bag

    return return_type(dsk, name, npartitions)


def _reduce(binop, sequence, initial=no_default):
    if initial is not no_default:
        return reduce(binop, sequence, initial)
    else:
        return reduce(binop, sequence)


def make_group(k, stage):
    def h(x):
        return x[0] // k ** stage % k
    return h


def groupby_tasks(b, grouper, hash=hash, max_branch=32):
    max_branch = max_branch or 32
    n = b.npartitions

    stages = int(math.ceil(math.log(n) / math.log(max_branch)))
    if stages > 1:
        k = int(math.ceil(n ** (1 / stages)))
    else:
        k = n

    groups = []
    splits = []
    joins = []

    inputs = [tuple(digit(i, j, k) for j in range(stages))
              for i in range(k**stages)]

    b2 = b.map(lambda x: (hash(grouper(x)), x))

    token = tokenize(b, grouper, hash, max_branch)

    start = dict((('shuffle-join-' + token, 0, inp),
                  (b2.name, i) if i < b.npartitions else [])
                 for i, inp in enumerate(inputs))

    for stage in range(1, stages + 1):
        group = dict((('shuffle-group-' + token, stage, inp),
                      (groupby,
                       (make_group, k, stage - 1),
                       ('shuffle-join-' + token, stage - 1, inp)))
                     for inp in inputs)

        split = dict((('shuffle-split-' + token, stage, i, inp),
                      (dict.get, ('shuffle-group-' + token, stage, inp), i, {}))
                     for i in range(k)
                     for inp in inputs)

        join = dict((('shuffle-join-' + token, stage, inp),
                     (list, (toolz.concat, [('shuffle-split-' + token, stage, inp[stage - 1],
                             insert(inp, stage - 1, j)) for j in range(k)])))
                    for inp in inputs)
        groups.append(group)
        splits.append(split)
        joins.append(join)

    end = dict((('shuffle-' + token, i),
                (list, (dict.items, (groupby, grouper, (pluck, 1, j)))))
               for i, j in enumerate(join))

    dsk = merge(b2.dask, start, end, *(groups + splits + joins))

    return type(b)(dsk, 'shuffle-' + token, len(inputs))


def groupby_disk(b, grouper, npartitions=None, blocksize=2**20):
    if npartitions is None:
        npartitions = b.npartitions
    token = tokenize(b, grouper, npartitions, blocksize)

    import partd
    p = ('partd-' + token,)
    dirname = _globals.get('temporary_directory', None)
    if dirname:
        file = (apply, partd.File, (), {'dir': dirname})
    else:
        file = (partd.File,)
    try:
        dsk1 = {p: (partd.Python, (partd.Snappy, file))}
    except AttributeError:
        dsk1 = {p: (partd.Python, file)}

    # Partition data on disk
    name = 'groupby-part-{0}-{1}'.format(funcname(grouper), token)
    dsk2 = dict(((name, i), (partition, grouper, (b.name, i),
                             npartitions, p, blocksize))
                for i in range(b.npartitions))

    # Barrier
    barrier_token = 'groupby-barrier-' + token

    def barrier(args):
        return 0

    dsk3 = {barrier_token: (barrier, list(dsk2))}

    # Collect groups
    name = 'groupby-collect-' + token
    dsk4 = dict(((name, i),
                 (collect, grouper, i, p, barrier_token))
                for i in range(npartitions))

    return type(b)(merge(b.dask, dsk1, dsk2, dsk3, dsk4), name, npartitions)


def empty_safe_apply(func, part, is_last):
    if isinstance(part, Iterator):
        try:
            _, part = peek(part)
        except StopIteration:
            if not is_last:
                return no_result
        return func(part)
    elif not is_last and len(part) == 0:
        return no_result
    else:
        return func(part)


def empty_safe_aggregate(func, parts, is_last):
    parts2 = (p for p in parts if p is not no_result)
    return empty_safe_apply(func, parts2, is_last)


def safe_take(n, b, warn=True):
    r = list(take(n, b))
    if len(r) != n and warn:
        warnings.warn("Insufficient elements for `take`. {0} elements "
                      "requested, only {1} elements available. Try passing "
                      "larger `npartitions` to `take`.".format(n, len(r)))
    return r


def random_sample(x, state_data, prob):
    """Filter elements of `x` by a probability `prob`.

    Parameters
    ----------
    x : iterable
    state_data : tuple
        A tuple that can be passed to ``random.Random``.
    prob : float
        A float between 0 and 1, representing the probability that each
        element will be yielded.
    """
    random_state = Random(state_data)
    for i in x:
        if random_state.random() < prob:
            yield i


def random_state_data_python(n, random_state=None):
    """Return a list of tuples that can initialize.
    ``random.Random``.

    Parameters
    ----------
    n : int
        Number of tuples to return.
    random_state : int or ``random.Random``, optional
        If an int, is used to seed a new ``random.Random``.
    """
    if not isinstance(random_state, Random):
        random_state = Random(random_state)

    maxuint32 = 1 << 32
    return [tuple(random_state.randint(0, maxuint32) for i in range(624))
            for i in range(n)]


def split(seq, n):
    """ Split apart a sequence into n equal pieces.

    >>> split(range(10), 3)
    [[0, 1, 2], [3, 4, 5], [6, 7, 8, 9]]
    """
    if not isinstance(seq, (list, tuple)):
        seq = list(seq)

    part = len(seq) / n
    L = [seq[int(part * i): int(part * (i + 1))] for i in range(n - 1)]
    L.append(seq[int(part * (n - 1)):])
    return L


def to_dataframe(seq, columns, dtypes):
    import pandas as pd
    seq = reify(seq)
    # pd.DataFrame expects lists, only copy if necessary
    if not isinstance(seq, list):
        seq = list(seq)
    res = pd.DataFrame(seq, columns=list(columns))
    return res.astype(dtypes, copy=False)
