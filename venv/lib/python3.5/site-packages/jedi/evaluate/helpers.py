import copy
import sys
import re
import os
from itertools import chain
from contextlib import contextmanager

from parso.python import tree

from jedi._compatibility import unicode
from jedi.parser_utils import get_parent_scope
from jedi.evaluate.compiled import CompiledObject


def is_stdlib_path(path):
    # Python standard library paths look like this:
    # /usr/lib/python3.5/...
    # TODO The implementation below is probably incorrect and not complete.
    if 'dist-packages' in path or 'site-packages' in path:
        return False

    base_path = os.path.join(sys.prefix, 'lib', 'python')
    return bool(re.match(re.escape(base_path) + '\d.\d', path))


def deep_ast_copy(obj):
    """
    Much, much faster than copy.deepcopy, but just for parser tree nodes.
    """
    # If it's already in the cache, just return it.
    new_obj = copy.copy(obj)

    # Copy children
    new_children = []
    for child in obj.children:
        if isinstance(child, tree.Leaf):
            new_child = copy.copy(child)
            new_child.parent = new_obj
        else:
            new_child = deep_ast_copy(child)
            new_child.parent = new_obj
        new_children.append(new_child)
    new_obj.children = new_children

    return new_obj


def evaluate_call_of_leaf(context, leaf, cut_own_trailer=False):
    """
    Creates a "call" node that consist of all ``trailer`` and ``power``
    objects.  E.g. if you call it with ``append``::

        list([]).append(3) or None

    You would get a node with the content ``list([]).append`` back.

    This generates a copy of the original ast node.

    If you're using the leaf, e.g. the bracket `)` it will return ``list([])``.

    We use this function for two purposes. Given an expression ``bar.foo``,
    we may want to
      - infer the type of ``foo`` to offer completions after foo
      - infer the type of ``bar`` to be able to jump to the definition of foo
    The option ``cut_own_trailer`` must be set to true for the second purpose.
    """
    trailer = leaf.parent
    # The leaf may not be the last or first child, because there exist three
    # different trailers: `( x )`, `[ x ]` and `.x`. In the first two examples
    # we should not match anything more than x.
    if trailer.type != 'trailer' or leaf not in (trailer.children[0], trailer.children[-1]):
        if trailer.type == 'atom':
            return context.eval_node(trailer)
        return context.eval_node(leaf)

    power = trailer.parent
    index = power.children.index(trailer)
    if cut_own_trailer:
        cut = index
    else:
        cut = index + 1

    if power.type == 'error_node':
        start = index
        while True:
            start -= 1
            base = power.children[start]
            if base.type != 'trailer':
                break
        trailers = power.children[start + 1: index + 1]
    else:
        base = power.children[0]
        trailers = power.children[1:cut]

    if base == 'await':
        base = trailers[0]
        trailers = trailers[1:]

    values = context.eval_node(base)
    from jedi.evaluate.syntax_tree import eval_trailer
    for trailer in trailers:
        values = eval_trailer(context, values, trailer)
    return values


def call_of_leaf(leaf):
    """
    Creates a "call" node that consist of all ``trailer`` and ``power``
    objects.  E.g. if you call it with ``append``::

        list([]).append(3) or None

    You would get a node with the content ``list([]).append`` back.

    This generates a copy of the original ast node.

    If you're using the leaf, e.g. the bracket `)` it will return ``list([])``.
    """
    # TODO this is the old version of this call. Try to remove it.
    trailer = leaf.parent
    # The leaf may not be the last or first child, because there exist three
    # different trailers: `( x )`, `[ x ]` and `.x`. In the first two examples
    # we should not match anything more than x.
    if trailer.type != 'trailer' or leaf not in (trailer.children[0], trailer.children[-1]):
        if trailer.type == 'atom':
            return trailer
        return leaf

    power = trailer.parent
    index = power.children.index(trailer)

    new_power = copy.copy(power)
    new_power.children = list(new_power.children)
    new_power.children[index + 1:] = []

    if power.type == 'error_node':
        start = index
        while True:
            start -= 1
            if power.children[start].type != 'trailer':
                break
        transformed = tree.Node('power', power.children[start:])
        transformed.parent = power.parent
        return transformed

    return power


def get_names_of_node(node):
    try:
        children = node.children
    except AttributeError:
        if node.type == 'name':
            return [node]
        else:
            return []
    else:
        return list(chain.from_iterable(get_names_of_node(c) for c in children))


def get_module_names(module, all_scopes):
    """
    Returns a dictionary with name parts as keys and their call paths as
    values.
    """
    names = chain.from_iterable(module.get_used_names().values())
    if not all_scopes:
        # We have to filter all the names that don't have the module as a
        # parent_scope. There's None as a parent, because nodes in the module
        # node have the parent module and not suite as all the others.
        # Therefore it's important to catch that case.
        names = [n for n in names if get_parent_scope(n).parent in (module, None)]
    return names


@contextmanager
def predefine_names(context, flow_scope, dct):
    predefined = context.predefined_names
    predefined[flow_scope] = dct
    try:
        yield
    finally:
        del predefined[flow_scope]


def is_compiled(context):
    return isinstance(context, CompiledObject)


def is_string(context):
    if context.evaluator.environment.version_info.major == 2:
        str_classes = (unicode, bytes)
    else:
        str_classes = (unicode,)
    return is_compiled(context) and isinstance(context.get_safe_value(default=None), str_classes)


def is_literal(context):
    return is_number(context) or is_string(context)


def _get_safe_value_or_none(context, accept):
    if is_compiled(context):
        value = context.get_safe_value(default=None)
        if isinstance(value, accept):
            return value


def get_int_or_none(context):
    return _get_safe_value_or_none(context, int)


def is_number(context):
    return _get_safe_value_or_none(context, (int, float)) is not None
