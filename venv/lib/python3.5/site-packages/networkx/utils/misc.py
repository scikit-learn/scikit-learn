"""
Miscellaneous Helpers for NetworkX.

These are not imported into the base networkx namespace but
can be accessed, for example, as

>>> import networkx
>>> networkx.utils.is_string_like('spam')
True
"""
# Authors:      Aric Hagberg (hagberg@lanl.gov),
#               Dan Schult(dschult@colgate.edu),
#               Ben Edwards(bedwards@cs.unm.edu)

#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
from collections import defaultdict
from collections import deque
import sys
import uuid
from itertools import tee, chain

# itertools.accumulate is only available on Python 3.2 or later.
#
# Once support for Python versions less than 3.2 is dropped, this code should
# be removed.
try:
    from itertools import accumulate
except ImportError:
    import operator

    # The code for this function is from the Python 3.5 documentation,
    # distributed under the PSF license:
    # <https://docs.python.org/3.5/library/itertools.html#itertools.accumulate>
    def accumulate(iterable, func=operator.add):
        it = iter(iterable)
        try:
            total = next(it)
        except StopIteration:
            return
        yield total
        for element in it:
            total = func(total, element)
            yield total

# 2.x/3.x compatibility
try:
    basestring
except NameError:
    basestring = str
    unicode = str

# some cookbook stuff
# used in deciding whether something is a bunch of nodes, edges, etc.
# see G.add_nodes and others in Graph Class in networkx/base.py


def is_string_like(obj):  # from John Hunter, types-free version
    """Check if obj is string."""
    return isinstance(obj, basestring)


def iterable(obj):
    """ Return True if obj is iterable with a well-defined len()."""
    if hasattr(obj, "__iter__"):
        return True
    try:
        len(obj)
    except:
        return False
    return True


def flatten(obj, result=None):
    """ Return flattened version of (possibly nested) iterable object. """
    if not iterable(obj) or is_string_like(obj):
        return obj
    if result is None:
        result = []
    for item in obj:
        if not iterable(item) or is_string_like(item):
            result.append(item)
        else:
            flatten(item, result)
    return obj.__class__(result)


def is_list_of_ints(intlist):
    """ Return True if list is a list of ints. """
    if not isinstance(intlist, list):
        return False
    for i in intlist:
        if not isinstance(i, int):
            return False
    return True


PY2 = sys.version_info[0] == 2
if PY2:
    def make_str(x):
        """Return the string representation of t."""
        if isinstance(x, unicode):
            return x
        else:
            # Note, this will not work unless x is ascii-encoded.
            # That is good, since we should be working with unicode anyway.
            # Essentially, unless we are reading a file, we demand that users
            # convert any encoded strings to unicode before using the library.
            #
            # Also, the str() is necessary to convert integers, etc.
            # unicode(3) works, but unicode(3, 'unicode-escape') wants a buffer.
            #
            return unicode(str(x), 'unicode-escape')
else:
    def make_str(x):
        """Return the string representation of t."""
        return str(x)


def generate_unique_node():
    """ Generate a unique node label."""
    return str(uuid.uuid1())


def default_opener(filename):
    """Opens `filename` using system's default program.

    Parameters
    ----------
    filename : str
        The path of the file to be opened.

    """
    from subprocess import call

    cmds = {'darwin': ['open'],
            'linux': ['xdg-open'],
            'linux2': ['xdg-open'],
            'win32': ['cmd.exe', '/C', 'start', '']}
    cmd = cmds[sys.platform] + [filename]
    call(cmd)


def dict_to_numpy_array(d, mapping=None):
    """Convert a dictionary of dictionaries to a numpy array
    with optional mapping."""
    try:
        return dict_to_numpy_array2(d, mapping)
    except (AttributeError, TypeError):
        # AttributeError is when no mapping was provided and v.keys() fails.
        # TypeError is when a mapping was provided and d[k1][k2] fails.
        return dict_to_numpy_array1(d, mapping)


def dict_to_numpy_array2(d, mapping=None):
    """Convert a dictionary of dictionaries to a 2d numpy array
    with optional mapping.

    """
    import numpy
    if mapping is None:
        s = set(d.keys())
        for k, v in d.items():
            s.update(v.keys())
        mapping = dict(zip(s, range(len(s))))
    n = len(mapping)
    a = numpy.zeros((n, n))
    for k1, i in mapping.items():
        for k2, j in mapping.items():
            try:
                a[i, j] = d[k1][k2]
            except KeyError:
                pass
    return a


def dict_to_numpy_array1(d, mapping=None):
    """Convert a dictionary of numbers to a 1d numpy array
    with optional mapping.

    """
    import numpy
    if mapping is None:
        s = set(d.keys())
        mapping = dict(zip(s, range(len(s))))
    n = len(mapping)
    a = numpy.zeros(n)
    for k1, i in mapping.items():
        i = mapping[k1]
        a[i] = d[k1]
    return a


def is_iterator(obj):
    """Returns True if and only if the given object is an iterator
    object.

    """
    has_next_attr = hasattr(obj, '__next__') or hasattr(obj, 'next')
    return iter(obj) is obj and has_next_attr


def arbitrary_element(iterable):
    """Returns an arbitrary element of `iterable` without removing it.

    This is most useful for "peeking" at an arbitrary element of a set,
    but can be used for any list, dictionary, etc., as well::

        >>> arbitrary_element({3, 2, 1})
        1
        >>> arbitrary_element('hello')
        'h'

    This function raises a :exc:`ValueError` if `iterable` is an
    iterator (because the current implementation of this function would
    consume an element from the iterator)::

        >>> iterator = iter([1, 2, 3])
        >>> arbitrary_element(iterator)
        Traceback (most recent call last):
            ...
        ValueError: cannot return an arbitrary item from an iterator

    """
    if is_iterator(iterable):
        raise ValueError('cannot return an arbitrary item from an iterator')
    # Another possible implementation is ``for x in iterable: return x``.
    return next(iter(iterable))


# Recipe from the itertools documentation.
def consume(iterator):
    "Consume the iterator entirely."
    # Feed the entire iterator into a zero-length deque.
    deque(iterator, maxlen=0)


# Recipe from the itertools documentation.
def pairwise(iterable, cyclic=False):
    "s -> (s0, s1), (s1, s2), (s2, s3), ..."
    a, b = tee(iterable)
    first = next(b, None)
    if cyclic is True:
        return zip(a, chain(b, (first,)))
    return zip(a, b)


def groups(many_to_one):
    """Converts a many-to-one mapping into a one-to-many mapping.

    `many_to_one` must be a dictionary whose keys and values are all
    :term:`hashable`.

    The return value is a dictionary mapping values from `many_to_one`
    to sets of keys from `many_to_one` that have that value.

    For example::

        >>> from networkx.utils import groups
        >>> many_to_one = {'a': 1, 'b': 1, 'c': 2, 'd': 3, 'e': 3}
        >>> groups(many_to_one)  # doctest: +SKIP
        {1: {'a', 'b'}, 2: {'c'}, 3: {'d', 'e'}}

    """
    one_to_many = defaultdict(set)
    for v, k in many_to_one.items():
        one_to_many[k].add(v)
    return dict(one_to_many)


def to_tuple(x):
    """Converts lists to tuples.

    For example::

        >>> from networkx.utils import to_tuple
        >>> a_list = [1, 2, [1, 4]]
        >>> to_tuple(a_list)
        (1, 2, (1, 4))

    """
    if not isinstance(x, (tuple, list)):
        return x
    return tuple(map(to_tuple, x))


def create_random_state(random_state=None):
    """Returns a numpy.random.RandomState instance depending on input.

    Parameters
    ----------
    random_state : int or RandomState instance or None  optional (default=None)
        If int, `random_state` is the seed used by the random number generator,
        if numpy.random.RandomState instance, `random_state` is the random
        number generator,
        if None, the random number generator is the RandomState instance used
        by numpy.random.
    """
    import numpy as np

    if random_state is None or random_state is np.random:
        return np.random.mtrand._rand
    if isinstance(random_state, np.random.RandomState):
        return random_state
    if isinstance(random_state, int):
        return np.random.RandomState(random_state)
    msg = '%r cannot be used to generate a numpy.random.RandomState instance'
    raise ValueError(msg % random_state)
