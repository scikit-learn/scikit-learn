"""
This module is not intended to be used in jedi, rather it will be fed to the
jedi-parser to replace classes in the typing module
"""

try:
    from collections import abc
except ImportError:
    # python 2
    import collections as abc


def factory(typing_name, indextypes):
    class Iterable(abc.Iterable):
        def __iter__(self):
            while True:
                yield indextypes[0]()

    class Iterator(Iterable, abc.Iterator):
        def next(self):
            """ needed for python 2 """
            return self.__next__()

        def __next__(self):
            return indextypes[0]()

    class Sequence(abc.Sequence):
        def __getitem__(self, index):
            return indextypes[0]()

    class MutableSequence(Sequence, abc.MutableSequence):
        pass

    class List(MutableSequence, list):
        pass

    class Tuple(Sequence, tuple):
        def __getitem__(self, index):
            if indextypes[1] == Ellipsis:
                # https://www.python.org/dev/peps/pep-0484/#the-typing-module
                # Tuple[int, ...] means a tuple of ints of indetermined length
                return indextypes[0]()
            else:
                return indextypes[index]()

    class AbstractSet(Iterable, abc.Set):
        pass

    class MutableSet(AbstractSet, abc.MutableSet):
        pass

    class KeysView(Iterable, abc.KeysView):
        pass

    class ValuesView(abc.ValuesView):
        def __iter__(self):
            while True:
                yield indextypes[1]()

    class ItemsView(abc.ItemsView):
        def __iter__(self):
            while True:
                yield (indextypes[0](), indextypes[1]())

    class Mapping(Iterable, abc.Mapping):
        def __getitem__(self, item):
            return indextypes[1]()

        def keys(self):
            return KeysView()

        def values(self):
            return ValuesView()

        def items(self):
            return ItemsView()

    class MutableMapping(Mapping, abc.MutableMapping):
        pass

    class Dict(MutableMapping, dict):
        pass

    dct = {
        "Sequence": Sequence,
        "MutableSequence": MutableSequence,
        "List": List,
        "Iterable": Iterable,
        "Iterator": Iterator,
        "AbstractSet": AbstractSet,
        "MutableSet": MutableSet,
        "Mapping": Mapping,
        "MutableMapping": MutableMapping,
        "Tuple": Tuple,
        "KeysView": KeysView,
        "ItemsView": ItemsView,
        "ValuesView": ValuesView,
        "Dict": Dict,
    }
    return dct[typing_name]
