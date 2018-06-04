from __future__ import absolute_import, division, print_function

import re
import os
from functools import partial

from .compatibility import apply
from .core import istask, get_dependencies, ishashable
from .utils import funcname, import_required


graphviz = import_required("graphviz", "Drawing dask graphs requires the "
                                       "`graphviz` python library and the "
                                       "`graphviz` system library to be "
                                       "installed.")


def task_label(task):
    """Label for a task on a dot graph.

    Examples
    --------
    >>> from operator import add
    >>> task_label((add, 1, 2))
    'add'
    >>> task_label((add, (add, 1, 2), 3))
    'add(...)'
    """
    func = task[0]
    if func is apply:
        func = task[1]
    if hasattr(func, 'funcs'):
        if len(func.funcs) > 1:
            return '{0}(...)'.format(funcname(func.funcs[0]))
        else:
            head = funcname(func.funcs[0])
    else:
        head = funcname(func)
    if any(has_sub_tasks(i) for i in task[1:]):
        return '{0}(...)'.format(head)
    else:
        return head


def has_sub_tasks(task):
    """Returns True if the task has sub tasks"""
    if istask(task):
        return True
    elif isinstance(task, list):
        return any(has_sub_tasks(i) for i in task)
    else:
        return False


def name(x):
    try:
        return str(hash(x))
    except TypeError:
        return str(hash(str(x)))


_HASHPAT = re.compile('([0-9a-z]{32})')
_UUIDPAT = re.compile('([0-9a-z]{8}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{12})')


def label(x, cache=None):
    """

    >>> label('x')
    'x'

    >>> label(('x', 1))
    "('x', 1)"

    >>> from hashlib import md5
    >>> x = 'x-%s-hello' % md5(b'1234').hexdigest()
    >>> x
    'x-81dc9bdb52d04dc20036dbd8313ed055-hello'

    >>> label(x)
    'x-#-hello'

    >>> from uuid import uuid1
    >>> x = 'x-%s-hello' % uuid1()
    >>> x  # doctest: +SKIP
    'x-4c1a3d7e-0b45-11e6-8334-54ee75105593-hello'

    >>> label(x)
    'x-#-hello'
    """
    s = str(x)
    for pattern in (_HASHPAT, _UUIDPAT):
        m = re.search(pattern, s)
        if m is not None:
            for h in m.groups():
                if cache is not None:
                    n = cache.get(h, len(cache))
                    label = '#{0}'.format(n)
                    # cache will be overwritten destructively
                    cache[h] = n
                else:
                    label = '#'
                s = s.replace(h, label)
    return s


def to_graphviz(dsk, data_attributes=None, function_attributes=None,
                rankdir='BT', graph_attr={}, node_attr=None, edge_attr=None, **kwargs):
    if data_attributes is None:
        data_attributes = {}
    if function_attributes is None:
        function_attributes = {}

    graph_attr = graph_attr or {}
    graph_attr['rankdir'] = rankdir
    graph_attr.update(kwargs)
    g = graphviz.Digraph(graph_attr=graph_attr,
                         node_attr=node_attr,
                         edge_attr=edge_attr)

    seen = set()
    cache = {}

    for k, v in dsk.items():
        k_name = name(k)
        if k_name not in seen:
            seen.add(k_name)
            attrs = data_attributes.get(k, {})
            attrs.setdefault('label', label(k, cache=cache))
            attrs.setdefault('shape', 'box')
            g.node(k_name, **attrs)

        if istask(v):
            func_name = name((k, 'function'))
            if func_name not in seen:
                seen.add(func_name)
                attrs = function_attributes.get(k, {})
                attrs.setdefault('label', task_label(v))
                attrs.setdefault('shape', 'circle')
                g.node(func_name, **attrs)
            g.edge(func_name, k_name)

            for dep in get_dependencies(dsk, k):
                dep_name = name(dep)
                if dep_name not in seen:
                    seen.add(dep_name)
                    attrs = data_attributes.get(dep, {})
                    attrs.setdefault('label', label(dep, cache=cache))
                    attrs.setdefault('shape', 'box')
                    g.node(dep_name, **attrs)
                g.edge(dep_name, func_name)
        elif ishashable(v) and v in dsk:
            g.edge(name(v), k_name)
    return g


IPYTHON_IMAGE_FORMATS = frozenset(['jpeg', 'png'])
IPYTHON_NO_DISPLAY_FORMATS = frozenset(['dot', 'pdf'])


def _get_display_cls(format):
    """
    Get the appropriate IPython display class for `format`.

    Returns `IPython.display.SVG` if format=='svg', otherwise
    `IPython.display.Image`.

    If IPython is not importable, return dummy function that swallows its
    arguments and returns None.
    """
    dummy = lambda *args, **kwargs: None
    try:
        import IPython.display as display
    except ImportError:
        # Can't return a display object if no IPython.
        return dummy

    if format in IPYTHON_NO_DISPLAY_FORMATS:
        # IPython can't display this format natively, so just return None.
        return dummy
    elif format in IPYTHON_IMAGE_FORMATS:
        # Partially apply `format` so that `Image` and `SVG` supply a uniform
        # interface to the caller.
        return partial(display.Image, format=format)
    elif format == 'svg':
        return display.SVG
    else:
        raise ValueError("Unknown format '%s' passed to `dot_graph`" % format)


def dot_graph(dsk, filename='mydask', format=None, **kwargs):
    """
    Render a task graph using dot.

    If `filename` is not None, write a file to disk with that name in the
    format specified by `format`.  `filename` should not include an extension.

    Parameters
    ----------
    dsk : dict
        The graph to display.
    filename : str or None, optional
        The name (without an extension) of the file to write to disk.  If
        `filename` is None, no file will be written, and we communicate with
        dot using only pipes.  Default is 'mydask'.
    format : {'png', 'pdf', 'dot', 'svg', 'jpeg', 'jpg'}, optional
        Format in which to write output file.  Default is 'png'.
    **kwargs
        Additional keyword arguments to forward to `to_graphviz`.

    Returns
    -------
    result : None or IPython.display.Image or IPython.display.SVG  (See below.)

    Notes
    -----
    If IPython is installed, we return an IPython.display object in the
    requested format.  If IPython is not installed, we just return None.

    We always return None if format is 'pdf' or 'dot', because IPython can't
    display these formats natively. Passing these formats with filename=None
    will not produce any useful output.

    See Also
    --------
    dask.dot.to_graphviz
    """
    g = to_graphviz(dsk, **kwargs)

    fmts = ['.png', '.pdf', '.dot', '.svg', '.jpeg', '.jpg']
    if format is None and any(filename.lower().endswith(fmt) for fmt in fmts):
        filename, format = os.path.splitext(filename)
        format = format[1:].lower()

    if format is None:
        format = 'png'

    data = g.pipe(format=format)
    if not data:
        raise RuntimeError("Graphviz failed to properly produce an image. "
                           "This probably means your installation of graphviz "
                           "is missing png support. See: "
                           "https://github.com/ContinuumIO/anaconda-issues/"
                           "issues/485 for more information.")

    display_cls = _get_display_cls(format)

    if not filename:
        return display_cls(data=data)

    full_filename = '.'.join([filename, format])
    with open(full_filename, 'wb') as f:
        f.write(data)

    return display_cls(filename=full_filename)
