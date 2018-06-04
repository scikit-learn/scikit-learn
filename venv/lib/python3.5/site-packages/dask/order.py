""" Static order of nodes in dask graph

Dask makes decisions on what tasks to prioritize both

*  Dynamically at runtime
*  Statically before runtime

Dynamically we prefer to run tasks that were just made available.  However when
several tasks become available at the same time we have an opportunity to break
ties in an intelligent way

        d
        |
    b   c
     \ /
      a

For example after we finish ``a`` we can choose to run either ``b`` or ``c``
next.  In this case we may prefer to start with ``c``, because it has other
dependents.

This is particularly important at the beginning of the computation when we
often dump hundreds of leaf nodes onto the scheduler at once.  The order in
which we start this computation can significantly change performance.


Breaking Ties
-------------

And so we create a total ordering over all nodes to serve as a tie breaker.  We
represent this ordering with a dictionary.  Lower scores have higher priority.

    {'d': 0,
     'c': 1,
     'a': 2,
     'b': 3}

There are several ways in which we might order our keys.  In practice we have
found the following objectives important:

1.  **Depth first**:  By traversing a tree deeply before broadly we encourage
    the completion of cohesive parts of the computation before starting new
    parts.  This helps to reduce memory footprint.
2.  **Lean**: We avoid tasks that have multiple dependents.  These tend to
    start new branches of computation that we would prefer to avoid until we
    have finished our current work.
3.  **Heavy-first**: When deciding between two possible dependencies we choose
    the dependency that is likely to take the longest first.  That way after it
    finishes it has to stay in memory only a short time while waiting for its
    co-dependent task to finish

So we perform a depth first search where we choose to traverse down children
with the following priority:

1.  Number of dependents (smaller is better)
2.  Number of total downstream dependencies (larger is better)
3.  Name of the task itself, as a tie breaker
"""
from __future__ import absolute_import, division, print_function

from .core import get_dependencies, reverse_dict, get_deps  # noqa: F401
from .utils_test import add, inc  # noqa: F401


def order(dsk, dependencies=None):
    """ Order nodes in dask graph

    This produces an ordering over our tasks that we use to break ties when
    executing.  We do this ahead of time to reduce a bit of stress on the
    scheduler and also to assist in static analysis.

    This currently traverses the graph as a single-threaded scheduler would
    traverse it.  It breaks ties in the following ways:

    1.  Start from roots nodes that have the largest subgraphs
    2.  When a node has dependencies that are not yet computed prefer
        dependencies with large subtrees  (start hard things first)
    2.  When we reach a node that can be computed we then traverse up and
        prefer dependents that have small super-trees (few total dependents)
        (finish existing work quickly)

    Examples
    --------
    >>> dsk = {'a': 1, 'b': 2, 'c': (inc, 'a'), 'd': (add, 'b', 'c')}
    >>> order(dsk)
    {'a': 0, 'c': 1, 'b': 2, 'd': 3}
    """
    if dependencies is None:
        dependencies = {k: get_dependencies(dsk, k) for k in dsk}

    for k, deps in dependencies.items():
        deps.discard(k)

    dependents = reverse_dict(dependencies)

    total_dependencies = ndependencies(dependencies, dependents)
    total_dependents = ndependents(dependencies, dependents)

    waiting = {k: set(v) for k, v in dependencies.items()}

    def dependencies_key(x):
        return total_dependencies.get(x, 0), ReverseStrComparable(x)

    def dependents_key(x):
        return (total_dependents.get(x, 0),
                StrComparable(x))

    result = dict()
    seen = set()  # tasks that should not be added again to the stack
    i = 0

    stack = [k for k, v in dependents.items() if not v]
    if len(stack) < 10000:
        stack = sorted(stack, key=dependencies_key)
    else:
        stack = stack[::-1]

    while stack:
        item = stack.pop()

        if item in result:
            continue

        deps = waiting[item]
        if deps:
            stack.append(item)
            seen.add(item)
            if len(deps) < 1000:
                deps = sorted(deps, key=dependencies_key)
            stack.extend(deps)
            continue

        result[item] = i
        i += 1

        for dep in dependents[item]:
            waiting[dep].discard(item)

        deps = [d for d in dependents[item]
                if d not in result and not (d in seen and waiting[d])]
        if len(deps) < 1000:
            deps = sorted(deps, key=dependents_key, reverse=True)
        stack.extend(deps)

    return result


def ndependents(dependencies, dependents):
    """ Number of total data elements that depend on key

    For each key we return the number of data that can only be run after this
    key is run.  The root nodes have value 1 while deep child nodes will have
    larger values.

    Examples
    --------
    >>> dsk = {'a': 1, 'b': (inc, 'a'), 'c': (inc, 'b')}
    >>> dependencies, dependents = get_deps(dsk)

    >>> sorted(ndependents(dependencies, dependents).items())
    [('a', 3), ('b', 2), ('c', 1)]
    """
    result = dict()
    num_needed = {k: len(v) for k, v in dependents.items()}
    current = {k for k, v in num_needed.items() if v == 0}
    while current:
        key = current.pop()
        result[key] = 1 + sum(result[parent] for parent in dependents[key])
        for child in dependencies[key]:
            num_needed[child] -= 1
            if num_needed[child] == 0:
                current.add(child)
    return result


def ndependencies(dependencies, dependents):
    """ Number of total data elements on which this key depends

    For each key we return the number of tasks that must be run for us to run
    this task.

    Examples
    --------
    >>> dsk = {'a': 1, 'b': (inc, 'a'), 'c': (inc, 'b')}
    >>> dependencies, dependents = get_deps(dsk)
    >>> sorted(ndependencies(dependencies, dependents).items())
    [('a', 1), ('b', 2), ('c', 3)]
    """
    result = dict()
    num_needed = {k: len(v) for k, v in dependencies.items()}
    current = {k for k, v in num_needed.items() if v == 0}
    while current:
        key = current.pop()
        result[key] = 1 + sum(result[child] for child in dependencies[key])
        for parent in dependents[key]:
            num_needed[parent] -= 1
            if num_needed[parent] == 0:
                current.add(parent)
    return result


class StrComparable(object):
    """ Wrap object so that it defaults to string comparison

    When comparing two objects of different types Python fails

    >>> 'a' < 1  # doctest: +SKIP
    Traceback (most recent call last):
        ...
    TypeError: '<' not supported between instances of 'str' and 'int'

    This class wraps the object so that, when this would occur it instead
    compares the string representation

    >>> StrComparable('a') < StrComparable(1)
    False
    """
    __slots__ = ('obj',)

    def __init__(self, obj):
        self.obj = obj

    def __lt__(self, other):
        try:
            return self.obj < other.obj
        except Exception:
            return str(self.obj) < str(other.obj)


class ReverseStrComparable(object):
    """ Wrap object so that it defaults to string comparison

    Used when sorting in reverse direction.  See StrComparable for normal
    documentation.
    """
    __slots__ = ('obj',)

    def __init__(self, obj):
        self.obj = obj

    def __lt__(self, other):
        try:
            return self.obj > other.obj
        except Exception:
            return str(self.obj) > str(other.obj)
