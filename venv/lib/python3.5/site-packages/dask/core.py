from __future__ import absolute_import, division, print_function

from itertools import chain

from .utils_test import add, inc  # noqa: F401


def ishashable(x):
    """ Is x hashable?

    Examples
    --------

    >>> ishashable(1)
    True
    >>> ishashable([1])
    False
    """
    try:
        hash(x)
        return True
    except TypeError:
        return False


def istask(x):
    """ Is x a runnable task?

    A task is a tuple with a callable first argument

    Examples
    --------

    >>> inc = lambda x: x + 1
    >>> istask((inc, 1))
    True
    >>> istask(1)
    False
    """
    return type(x) is tuple and x and callable(x[0])


def has_tasks(dsk, x):
    """Whether ``x`` has anything to compute.

    Returns True if:
    - ``x`` is a task
    - ``x`` is a key in ``dsk``
    - ``x`` is a list that contains any tasks or keys
    """
    if istask(x):
        return True
    try:
        if x in dsk:
            return True
    except Exception:
        pass
    if isinstance(x, list):
        for i in x:
            if has_tasks(dsk, i):
                return True
    return False


def preorder_traversal(task):
    """A generator to preorder-traverse a task."""

    for item in task:
        if istask(item):
            for i in preorder_traversal(item):
                yield i
        elif isinstance(item, list):
            yield list
            for i in preorder_traversal(item):
                yield i
        else:
            yield item


def _get_nonrecursive(d, x, maxdepth=1000):
    # Non-recursive. DAG property is checked upon reaching maxdepth.
    _list = lambda *args: list(args)

    # We construct a nested hierarchy of tuples to mimic the execution stack
    # of frames that Python would maintain for a recursive implementation.
    # A frame is associated with a single task from a Dask.
    # A frame tuple has three elements:
    #    1) The function for the task.
    #    2) The arguments for the task (typically keys in the Dask).
    #       Arguments are stored in reverse order, and elements are popped
    #       as they are evaluated.
    #    3) The calculated results of the arguments from (2).
    stack = [(lambda x: x, [x], [])]
    while True:
        func, args, results = stack[-1]
        if not args:
            val = func(*results)
            if len(stack) == 1:
                return val
            stack.pop()
            stack[-1][2].append(val)
            continue
        elif maxdepth and len(stack) > maxdepth:
            cycle = getcycle(d, x)
            if cycle:
                cycle = '->'.join(cycle)
                raise RuntimeError('Cycle detected in Dask: %s' % cycle)
            maxdepth = None

        key = args.pop()
        if isinstance(key, list):
            stack.append((_list, list(key[::-1]), []))
            continue
        elif ishashable(key) and key in d:
            args.append(d[key])
            continue
        elif istask(key):
            stack.append((key[0], list(key[:0:-1]), []))
        else:
            results.append(key)


def _get_recursive(d, x):
    # recursive, no cycle detection
    if isinstance(x, list):
        return [_get_recursive(d, k) for k in x]
    elif ishashable(x) and x in d:
        return _get_recursive(d, d[x])
    elif istask(x):
        func, args = x[0], x[1:]
        args2 = [_get_recursive(d, k) for k in args]
        return func(*args2)
    else:
        return x


def get(d, x, recursive=False):
    """ Get value from Dask

    Examples
    --------

    >>> inc = lambda x: x + 1
    >>> d = {'x': 1, 'y': (inc, 'x')}

    >>> get(d, 'x')
    1
    >>> get(d, 'y')
    2
    """
    _get = _get_recursive if recursive else _get_nonrecursive
    if isinstance(x, list):
        return tuple(get(d, k) for k in x)
    elif x in d:
        return _get(d, x)
    raise KeyError("{0} is not a key in the graph".format(x))


def get_dependencies(dsk, key=None, task=None, as_list=False):
    """ Get the immediate tasks on which this task depends

    Examples
    --------
    >>> dsk = {'x': 1,
    ...        'y': (inc, 'x'),
    ...        'z': (add, 'x', 'y'),
    ...        'w': (inc, 'z'),
    ...        'a': (add, (inc, 'x'), 1)}

    >>> get_dependencies(dsk, 'x')
    set([])

    >>> get_dependencies(dsk, 'y')
    set(['x'])

    >>> get_dependencies(dsk, 'z')  # doctest: +SKIP
    set(['x', 'y'])

    >>> get_dependencies(dsk, 'w')  # Only direct dependencies
    set(['z'])

    >>> get_dependencies(dsk, 'a')  # Ignore non-keys
    set(['x'])

    >>> get_dependencies(dsk, task=(inc, 'x'))  # provide tasks directly
    set(['x'])
    """
    if key is not None:
        arg = dsk[key]
    elif task is not None:
        arg = task
    else:
        raise ValueError("Provide either key or task")

    result = []
    work = [arg]

    while work:
        new_work = []
        for w in work:
            typ = type(w)
            if typ is tuple and w and callable(w[0]):  # istask(w)
                new_work += w[1:]
            elif typ is list:
                new_work += w
            elif typ is dict:
                new_work += w.values()
            else:
                try:
                    if w in dsk:
                        result.append(w)
                except TypeError:  # not hashable
                    pass
        work = new_work

    return result if as_list else set(result)


def get_deps(dsk):
    """ Get dependencies and dependents from dask dask graph

    >>> dsk = {'a': 1, 'b': (inc, 'a'), 'c': (inc, 'b')}
    >>> dependencies, dependents = get_deps(dsk)
    >>> dependencies
    {'a': set([]), 'c': set(['b']), 'b': set(['a'])}
    >>> dependents
    {'a': set(['b']), 'c': set([]), 'b': set(['c'])}
    """
    dependencies = {k: get_dependencies(dsk, task=v)
                    for k, v in dsk.items()}
    dependents = reverse_dict(dependencies)
    return dependencies, dependents


def flatten(seq, container=list):
    """

    >>> list(flatten([1]))
    [1]

    >>> list(flatten([[1, 2], [1, 2]]))
    [1, 2, 1, 2]

    >>> list(flatten([[[1], [2]], [[1], [2]]]))
    [1, 2, 1, 2]

    >>> list(flatten(((1, 2), (1, 2)))) # Don't flatten tuples
    [(1, 2), (1, 2)]

    >>> list(flatten((1, 2, [3, 4]))) # support heterogeneous
    [1, 2, 3, 4]
    """
    if isinstance(seq, str):
        yield seq
    else:
        for item in seq:
            if isinstance(item, container):
                for item2 in flatten(item, container=container):
                    yield item2
            else:
                yield item


def reverse_dict(d):
    """

    >>> a, b, c = 'abc'
    >>> d = {a: [b, c], b: [c]}
    >>> reverse_dict(d)  # doctest: +SKIP
    {'a': set([]), 'b': set(['a']}, 'c': set(['a', 'b'])}
    """
    terms = list(d.keys()) + list(chain.from_iterable(d.values()))
    result = {t: set() for t in terms}
    for k, vals in d.items():
        for val in vals:
            result[val].add(k)
    return result


def subs(task, key, val):
    """ Perform a substitution on a task

    Examples
    --------

    >>> subs((inc, 'x'), 'x', 1)  # doctest: +SKIP
    (inc, 1)
    """
    type_task = type(task)
    if not (type_task is tuple and task and callable(task[0])):  # istask(task):
        try:
            if type_task is type(key) and task == key:
                return val
        except Exception:
            pass
        if type_task is list:
            return [subs(x, key, val) for x in task]
        return task
    newargs = []
    for arg in task[1:]:
        type_arg = type(arg)
        if type_arg is tuple and arg and callable(arg[0]):  # istask(task):
            arg = subs(arg, key, val)
        elif type_arg is list:
            arg = [subs(x, key, val) for x in arg]
        elif type_arg is type(key):
            try:
                # Can't do a simple equality check, since this may trigger
                # a FutureWarning from NumPy about array equality
                # https://github.com/dask/dask/pull/2457
                if len(arg) == len(key) and all(type(aa) == type(bb) and aa == bb
                                                for aa, bb in zip(arg, key)):
                    arg = val

            except (TypeError, AttributeError):
                # Handle keys which are not sized (len() fails), but are hashable
                if arg == key:
                    arg = val
        newargs.append(arg)
    return task[:1] + tuple(newargs)


def _toposort(dsk, keys=None, returncycle=False, dependencies=None):
    # Stack-based depth-first search traversal.  This is based on Tarjan's
    # method for topological sorting (see wikipedia for pseudocode)
    if keys is None:
        keys = dsk
    elif not isinstance(keys, list):
        keys = [keys]
    if not returncycle:
        ordered = []

    # Nodes whose descendents have been completely explored.
    # These nodes are guaranteed to not be part of a cycle.
    completed = set()

    # All nodes that have been visited in the current traversal.  Because
    # we are doing depth-first search, going "deeper" should never result
    # in visiting a node that has already been seen.  The `seen` and
    # `completed` sets are mutually exclusive; it is okay to visit a node
    # that has already been added to `completed`.
    seen = set()

    if dependencies is None:
        dependencies = dict((k, get_dependencies(dsk, k)) for k in dsk)

    for key in keys:
        if key in completed:
            continue
        nodes = [key]
        while nodes:
            # Keep current node on the stack until all descendants are visited
            cur = nodes[-1]
            if cur in completed:
                # Already fully traversed descendants of cur
                nodes.pop()
                continue
            seen.add(cur)

            # Add direct descendants of cur to nodes stack
            next_nodes = []
            for nxt in dependencies[cur]:
                if nxt not in completed:
                    if nxt in seen:
                        # Cycle detected!
                        cycle = [nxt]
                        while nodes[-1] != nxt:
                            cycle.append(nodes.pop())
                        cycle.append(nodes.pop())
                        cycle.reverse()
                        if returncycle:
                            return cycle
                        else:
                            cycle = '->'.join(str(x) for x in cycle)
                            raise RuntimeError('Cycle detected in Dask: %s' % cycle)
                    next_nodes.append(nxt)

            if next_nodes:
                nodes.extend(next_nodes)
            else:
                # cur has no more descendants to explore, so we're done with it
                if not returncycle:
                    ordered.append(cur)
                completed.add(cur)
                seen.remove(cur)
                nodes.pop()
    if returncycle:
        return []
    return ordered


def toposort(dsk, dependencies=None):
    """ Return a list of keys of dask sorted in topological order."""
    return _toposort(dsk, dependencies=dependencies)


def getcycle(d, keys):
    """ Return a list of nodes that form a cycle if Dask is not a DAG.

    Returns an empty list if no cycle is found.

    ``keys`` may be a single key or list of keys.

    Examples
    --------

    >>> d = {'x': (inc, 'z'), 'y': (inc, 'x'), 'z': (inc, 'y')}
    >>> getcycle(d, 'x')
    ['x', 'z', 'y', 'x']

    See Also
    --------
    isdag
    """
    return _toposort(d, keys=keys, returncycle=True)


def isdag(d, keys):
    """ Does Dask form a directed acyclic graph when calculating keys?

    ``keys`` may be a single key or list of keys.

    Examples
    --------

    >>> inc = lambda x: x + 1
    >>> isdag({'x': 0, 'y': (inc, 'x')}, 'y')
    True
    >>> isdag({'x': (inc, 'y'), 'y': (inc, 'x')}, 'y')
    False

    See Also
    --------
    getcycle
    """
    return not getcycle(d, keys)


class literal(object):
    """A small serializable object to wrap literal values without copying"""
    __slots__ = ('data',)

    def __init__(self, data):
        self.data = data

    def __repr__(self):
        return 'literal<type=%s>' % type(self.data).__name__

    def __reduce__(self):
        return (literal, (self.data,))

    def __call__(self):
        return self.data


def quote(x):
    """ Ensure that this value remains this value in a dask graph

    Some values in dask graph take on special meaning. Sometimes we want to
    ensure that our data is not interpreted but remains literal.

    >>> quote((add, 1, 2))  # doctest: +SKIP
    (literal<type=tuple>,)
    """
    if istask(x) or type(x) is list:
        return (literal(x),)
    return x
