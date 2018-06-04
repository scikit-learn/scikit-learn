from toolz import concat, unique, count
from collections import Mapping


class ShareDict(Mapping):
    """ A Mapping composed of other Mappings

    This is a union of other disjoint mappings.  It allows the combination of
    many dicts into a single dict-like object without creating copies of the
    underlying dicts.  It provides cheap ``update``, ``len`` and ``__iter__``
    operations as well as a fairly cheap ``__getitem__`` operation (linear in
    the number of constituent mappings).

    This class is optimized for Dask's use, and may not be generally useful.
    Users may want to consider the standard ``collections.ChainMap`` data
    structure.

    This class makes the following assumptions:

    1.  Constituent mappings are disjoint.  No key is in more than one mapping.
    2.  Constituent mappings will not be modified

    Note that ShareDict does not enforce these assumptions.  It is up to the
    user to guarantee them.

    Examples
    --------
    >>> a = {'x': 1, 'y': 2}
    >>> b = {'z': 3}
    >>> s = ShareDict()
    >>> s.update(a)
    >>> s.update(b)

    >>> dict(s)  # doctest: +SKIP
    {'x': 1, 'y': 2, 'z': 3}

    These dictionaries are stored within an internal dictionary of dictionaries

    >>> list(s.dicts.values())  # doctest: +SKIP
    [{'x': 1, 'y': 2}, {'z': 3}]

    By default these are named by their object id.  However, you can also
    provide explicit names.

    >>> s = ShareDict()
    >>> s.update_with_key(a, key='a')
    >>> s.update_with_key(b, key='b')
    >>> s.dicts  # doctest: +SKIP
    {'a': {'x': 1, 'y': 2}, 'b': {'z': 3}}
    """
    def __init__(self):
        self.dicts = dict()

    def update_with_key(self, arg, key=None):
        if type(arg) is ShareDict:
            assert key is None
            self.dicts.update(arg.dicts)
            return

        if key is None:
            key = id(arg)

        assert isinstance(arg, dict)
        if arg:
            self.dicts[key] = arg

    def update(self, arg):
        self.update_with_key(arg)

    def __getitem__(self, key):
        for d in self.dicts.values():
            if key in d:
                return d[key]
        raise KeyError(key)

    def __len__(self):
        return count(iter(self))

    def items(self):
        seen = set()
        for d in self.dicts.values():
            for key in d:
                if key not in seen:
                    seen.add(key)
                    yield (key, d[key])

    def __iter__(self):
        return unique(concat(self.dicts.values()))


def merge(*dicts):
    result = ShareDict()
    for d in dicts:
        if isinstance(d, tuple):
            key, d = d
            result.update_with_key(d, key=key)
        else:
            result.update_with_key(d)
    return result
